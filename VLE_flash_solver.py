import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union
from scipy.optimize import brentq


# ==============================================================================
# STEP 1A: PURE-COMPONENT VAPOR PRESSURE LOADER (PerryData.xlsx + Antoine)
# ==============================================================================

@dataclass
class PureComponent:
    """
    Generalized pure-component thermodynamic container.
    Supports either:
      - DIPPR 101 coefficients from PerryData.xlsx:
            ln(P_Pa) = C1 + C2/T + C3*ln(T) + C4*(T**C5)   [T in K, P in Pa]
      - Antoine coefficients (optional or standalone):
            log10(P_mmHg) = A - B / (C + T_C)             [T in °C, P in mmHg]
    """
    name: str
    MW: float
    dippr_coeffs: Optional[Tuple[float, float, float, float, float]] = None
    Tmin_K: float = 150.0
    Tmax_K: float = 650.0
    antoine_ABC: Optional[Tuple[float, float, float]] = None

    def vapor_pressure_Pa(self, T_C: float, enforce_bounds: bool = True) -> float:
        """Compute pure-component saturation vapor pressure in Pascals (Pa)."""
        T_K = T_C + 273.15
        if self.dippr_coeffs is not None:
            if enforce_bounds and not (self.Tmin_K <= T_K <= self.Tmax_K):
                raise ValueError(
                    f"{self.name}: T={T_K:.2f} K is outside PerryData valid range "
                    f"[{self.Tmin_K:.2f}, {self.Tmax_K:.2f}] K."
                )
            C1, C2, C3, C4, C5 = self.dippr_coeffs
            return float(np.exp(C1 + C2 / T_K + C3 * np.log(T_K) + C4 * (T_K ** C5)))

        if self.antoine_ABC is not None:
            if enforce_bounds and not (-50.0 <= T_C <= 350.0):
                raise ValueError(f"{self.name}: T={T_C:.2f} °C is outside Antoine range.")
            A, B, C = self.antoine_ABC
            P_mmHg = 10.0 ** (A - B / (C + T_C))
            return float(P_mmHg * 133.322)

        raise ValueError(f"No vapor pressure correlation available for '{self.name}'.")

    def saturation_temperature_C(self, P_Pa: float) -> float:
        """Compute normal/saturation boiling point (°C) at pressure P_Pa."""
        if self.dippr_coeffs is not None:
            return float(brentq(
                lambda t: self.vapor_pressure_Pa(t, enforce_bounds=False) - P_Pa,
                self.Tmin_K - 273.15,
                self.Tmax_K - 273.15,
                xtol=1e-8,
            ))
        if self.antoine_ABC is not None:
            A, B, C = self.antoine_ABC
            P_mmHg = P_Pa / 133.322
            return float(B / (A - np.log10(P_mmHg)) - C)

        raise ValueError(f"No vapor pressure correlation available for '{self.name}'.")


def load_pure_components(
    requested_names: Dict[str, str],
    perry_excel_path: Optional[str] = None,
    custom_components: Optional[Dict[str, PureComponent]] = None,
) -> Dict[str, PureComponent]:
    """
    ============================================================================
    [EXCEL IMPORT BOUNDARY — PerryData.xlsx (Optional)]
    Builds a PureComponent dictionary from:
      1) `custom_components` (user-supplied Antoine or DIPPR PureComponent objects), and/or
      2) `perry_excel_path` ('Pvap' + 'LiqDens' sheets from PerryData.xlsx).
    Users can supply ONLY Antoine constants by leaving `perry_excel_path=None`.
    ============================================================================
    """
    db: Dict[str, PureComponent] = {}
    merged = None

    if perry_excel_path is not None:
        # <-- EXCEL IMPORT: Reading 'Pvap' and 'LiqDens' sheets from PerryData.xlsx
        pvap_df = pd.read_excel(perry_excel_path, sheet_name="Pvap")
        dens_df = pd.read_excel(perry_excel_path, sheet_name="LiqDens")[["CmpdNo", "MolWt"]]
        merged = pd.merge(pvap_df, dens_df, on="CmpdNo", how="left")

    for key, query_name in requested_names.items():
        if custom_components and key in custom_components:
            db[key] = custom_components[key]
            continue

        if merged is None:
            raise ValueError(
                f"Component '{key}' not found in `custom_components`, and no `perry_excel_path` was provided."
            )

        mask = (
            merged["Name"].astype(str).str.strip().str.lower() == str(query_name).strip().lower()
        ) | (
            merged["CASNo"].astype(str).str.strip() == str(query_name).strip()
        )
        matches = merged[mask]
        if matches.empty:
            raise KeyError(
                f"Component '{query_name}' (key='{key}') not found in '{perry_excel_path}' (Pvap sheet). "
                f"Supply it via `custom_components`."
            )
        row = matches.iloc[0]
        mw_val = float(row["MolWt"]) if pd.notna(row["MolWt"]) else float("nan")
        db[key] = PureComponent(
            name=str(row["Name"]).strip(),
            MW=mw_val,
            dippr_coeffs=(
                float(row["C1"]),
                float(row["C2"]),
                float(row["C3"]),
                float(row["C4"]),
                float(row["C5"]),
            ),
            Tmin_K=float(row["Tmin"]),
            Tmax_K=float(row["Tmax K"]),
        )
    return db


# ==============================================================================
# STEP 1B: GENERALIZED MODIFIED UNIFAC (DORTMUND) DATABASE LOADER & MODEL
# ==============================================================================

SubgroupID = Union[str, Tuple[int, str]]  # e.g. "ACH" or (13, "CHO")


class ModUNIFACDatabase:
    """
    ============================================================================
    [EXCEL IMPORT BOUNDARY — MODUNI.xlsx]
    Imports Modified UNIFAC (Dortmund) subgroup volume/surface parameters (R, Q)
    from the 'Groups' sheet and main-group interaction parameter matrices from
    sheets 'A', 'B', and 'C'. See README.md for exact sheet layouts and units.
    ============================================================================
    """

    def __init__(self, moduni_excel_path: str):
        self.excel_path = moduni_excel_path
        self._load_excel_sheets()

    def _load_excel_sheets(self):
        # <-- EXCEL IMPORT: Reading subgroup R_k and Q_k table from 'Groups' sheet
        raw_groups = pd.read_excel(self.excel_path, sheet_name="Groups", skiprows=1)
        raw_groups = raw_groups.drop(
            columns=[c for c in raw_groups.columns if "Unnamed" in str(c)]
        )
        raw_groups.columns = ["sr_no", "subgrp", "main_grp_no", "main_grp_name", "R", "Q"]
        raw_groups = raw_groups.dropna(subset=["subgrp"]).copy()

        raw_groups["sr_no"] = raw_groups["sr_no"].astype(int)
        raw_groups["subgrp"] = raw_groups["subgrp"].astype(str).str.strip()
        raw_groups["main_grp_no"] = raw_groups["main_grp_no"].astype(int)
        raw_groups["R"] = raw_groups["R"].astype(float)
        raw_groups["Q"] = raw_groups["Q"].astype(float)
        self.groups_df = raw_groups

        self._by_pair: Dict[Tuple[int, str], dict] = {}
        self._by_name: Dict[str, List[dict]] = {}

        for _, r in raw_groups.iterrows():
            entry = {
                "sr_no": int(r["sr_no"]),
                "subgrp": r["subgrp"],
                "main_grp_no": int(r["main_grp_no"]),
                "main_grp_name": str(r["main_grp_name"]).strip(),
                "R": float(r["R"]),
                "Q": float(r["Q"]),
            }
            self._by_pair[(entry["main_grp_no"], entry["subgrp"])] = entry
            self._by_name.setdefault(entry["subgrp"], []).append(entry)

        # <-- EXCEL IMPORT: Reading main-group interaction matrices A, B, C
        self.A = self._read_matrix_sheet("A")
        self.B = self._read_matrix_sheet("B")
        self.C = self._read_matrix_sheet("C")

    def _read_matrix_sheet(self, sheet_name: str) -> pd.DataFrame:
        # <-- EXCEL IMPORT: Matrix sheet ('A', 'B', or 'C') indexed by Main Group ID
        raw = pd.read_excel(self.excel_path, sheet_name=sheet_name, header=None)
        col_ids = raw.iloc[1, 2:].astype(int).values
        row_ids = raw.iloc[2:, 1].astype(int).values
        vals = raw.iloc[2:, 2:].astype(float).values
        return pd.DataFrame(vals, index=row_ids, columns=col_ids)

    def resolve_subgroup(self, sg_id: SubgroupID) -> dict:
        """Resolve a subgroup identifier ('ACH' or (13, 'CHO')) to its parameter dict."""
        if isinstance(sg_id, tuple):
            if sg_id not in self._by_pair:
                raise KeyError(
                    f"Subgroup (main_grp={sg_id[0]}, subgrp='{sg_id[1]}') not found in MODUNI.xlsx."
                )
            return self._by_pair[sg_id]

        matches = self._by_name.get(sg_id.strip(), [])
        if not matches:
            raise KeyError(f"Subgroup '{sg_id}' not found in MODUNI.xlsx.")
        if len(matches) > 1:
            options = [(m["main_grp_no"], m["subgrp"], m["main_grp_name"]) for m in matches]
            raise ValueError(
                f"Subgroup '{sg_id}' is ambiguous in MODUNI.xlsx (matches {options}). "
                f"Specify it as a tuple (main_grp_no, '{sg_id}'), e.g., {(matches[0]['main_grp_no'], sg_id)}."
            )
        return matches[0]


class ModUNIFACModel:
    """
    Generalized multicomponent Modified UNIFAC (Dortmund) activity coefficient model
    (Weidlich & Gmehling, 1987; Gmehling et al., 1993).
    """

    def __init__(
        self,
        component_keys: List[str],
        component_subgroups: Dict[str, Dict[SubgroupID, float]],
        db: ModUNIFACDatabase,
        T_K: float,
    ):
        self.keys = component_keys
        self.n_comp = len(component_keys)
        self.db = db
        self.T_K = float(T_K)

        unique_sg_entries: List[dict] = []
        sg_index_map: Dict[Tuple[int, str], int] = {}

        for ck in component_keys:
            if ck not in component_subgroups:
                raise KeyError(f"Subgroup decomposition missing for component '{ck}'.")
            for sg_id in component_subgroups[ck]:
                entry = db.resolve_subgroup(sg_id)
                pair_key = (entry["main_grp_no"], entry["subgrp"])
                if pair_key not in sg_index_map:
                    sg_index_map[pair_key] = len(unique_sg_entries)
                    unique_sg_entries.append(entry)

        self.n_sub = len(unique_sg_entries)
        self.R_k = np.array([e["R"] for e in unique_sg_entries], dtype=float)
        self.Q_k = np.array([e["Q"] for e in unique_sg_entries], dtype=float)
        self.main_grps = [e["main_grp_no"] for e in unique_sg_entries]

        # Stoichiometric subgroup matrix nu[i, k]
        self.nu = np.zeros((self.n_comp, self.n_sub), dtype=float)
        for i, ck in enumerate(component_keys):
            for sg_id, count in component_subgroups[ck].items():
                entry = db.resolve_subgroup(sg_id)
                k = sg_index_map[(entry["main_grp_no"], entry["subgrp"])]
                self.nu[i, k] = float(count)

        self.r = self.nu @ self.R_k
        self.q = self.nu @ self.Q_k
        self._build_temperature_terms()

    def _build_temperature_terms(self):
        T = self.T_K
        self.Psi = np.zeros((self.n_sub, self.n_sub), dtype=float)

        for n in range(self.n_sub):
            gn = self.main_grps[n]
            for m in range(self.n_sub):
                gm = self.main_grps[m]
                a_nm = self.db.A.loc[gn, gm]
                b_nm = self.db.B.loc[gn, gm]
                c_nm = self.db.C.loc[gn, gm]
                if pd.isna(a_nm) or pd.isna(b_nm) or pd.isna(c_nm):
                    raise ValueError(
                        f"Missing Modified UNIFAC interaction parameter in MODUNI.xlsx "
                        f"between Main Group {gn} and Main Group {gm}."
                    )
                self.Psi[n, m] = np.exp(-(a_nm + b_nm * T + c_nm * (T ** 2)) / T)

        self.ln_Gamma_pure = np.zeros((self.n_comp, self.n_sub), dtype=float)
        for i in range(self.n_comp):
            x_pure = np.zeros(self.n_comp, dtype=float)
            x_pure[i] = 1.0
            self.ln_Gamma_pure[i] = self._group_residual_ln_gamma(x_pure)

    def _group_residual_ln_gamma(self, x: np.ndarray) -> np.ndarray:
        group_moles = x @ self.nu
        X_k = group_moles / np.sum(group_moles)
        Theta = (self.Q_k * X_k) / np.sum(self.Q_k * X_k)
        S = Theta @ self.Psi
        term = self.Psi @ (Theta / S)
        return self.Q_k * (1.0 - np.log(S) - term)

    def activity_coefficients(self, x: np.ndarray) -> np.ndarray:
        """Compute activity coefficients gamma_i for all components at composition x."""
        x = np.maximum(np.asarray(x, dtype=float), 1e-15)
        x = x / np.sum(x)

        # 1. Combinatorial term (Weidlich & Gmehling 3/4-exponent modification)
        r_34 = self.r ** 0.75
        V_prime = r_34 / np.sum(x * r_34)
        V = self.r / np.sum(x * self.r)
        F = self.q / np.sum(x * self.q)
        ln_gamma_C = 1.0 - V_prime + np.log(V_prime) - 5.0 * self.q * (
            1.0 - (V / F) + np.log(V / F)
        )

        # 2. Residual term
        ln_Gamma_mix = self._group_residual_ln_gamma(x)
        ln_gamma_R = np.sum(self.nu * (ln_Gamma_mix - self.ln_Gamma_pure), axis=1)

        return np.exp(ln_gamma_C + ln_gamma_R)


# ==============================================================================
# STEP 2A: RACHFORD-RICE SOLVER
# ==============================================================================

def rachford_rice(z: np.ndarray, K: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Solve the Rachford-Rice equation for overall vapour fraction beta:
        sum_i [ z_i * (K_i - 1) / (1 + beta * (K_i - 1)) ] = 0
    """
    if np.all(K <= 1.0):
        return 0.0, z.copy(), K * z / (K * z).sum() if (K * z).sum() > 0 else z.copy()
    if np.all(K >= 1.0):
        x_trivial = z / K
        x_trivial /= x_trivial.sum()
        return 1.0, x_trivial, z.copy()

    def rr_func(beta):
        return np.sum(z * (K - 1.0) / (1.0 + beta * (K - 1.0)))

    # Michelsen (1982) bounds: 1/(1 - K_max) < beta < 1/(1 - K_min)
    beta_lo = max(0.0, 1.0 / (1.0 - np.max(K)) + 1e-10)
    beta_hi = min(1.0, 1.0 / (1.0 - np.min(K)) - 1e-10)

    try:
        beta = brentq(rr_func, beta_lo, beta_hi, xtol=1e-12, maxiter=200)
    except ValueError:
        beta = 0.5

    denom = 1.0 + beta * (K - 1.0)
    x = z / denom
    x /= x.sum()
    y = K * x
    y /= y.sum()
    return float(beta), x, y


# ==============================================================================
# STEP 2B: GENERALIZED MODIFIED UNIFAC FLASH SOLVER
# ==============================================================================

class FlashSolver:
    """
    Isothermal-isobaric (PT) multicomponent VLE flash calculator powered
    exclusively by Modified UNIFAC (Dortmund) activity coefficients.
    """

    def __init__(
        self,
        component_keys: List[str],
        component_db: Dict[str, PureComponent],
        component_subgroups: Dict[str, Dict[SubgroupID, float]],
        moduni_db: ModUNIFACDatabase,
        T_C: float,
        P_Pa: float = 101325.0,
    ):
        for k in component_keys:
            if k not in component_db:
                raise KeyError(
                    f"Component '{k}' not found in `component_db`. "
                    f"Available: {list(component_db.keys())}"
                )
        self.keys = component_keys
        self.comps = [component_db[k] for k in component_keys]
        self.subgroups = component_subgroups
        self.moduni_db = moduni_db
        self.n = len(component_keys)
        self.T_C = float(T_C)
        self.T_K = self.T_C + 273.15
        self.P_Pa = float(P_Pa)
        self._unifac = ModUNIFACModel(component_keys, component_subgroups, moduni_db, self.T_K)

    # ── State Setters ────────────────────────────────────────────────────────

    def set_pressure(self, P_Pa: float):
        """Update operating pressure (Pa) without rebuilding UNIFAC matrices."""
        self.P_Pa = float(P_Pa)

    def set_temperature(self, T_C: float):
        """Update operating temperature (°C) and rebuild temperature-dependent Psi_nm."""
        self.T_C = float(T_C)
        self.T_K = self.T_C + 273.15
        self._unifac = ModUNIFACModel(self.keys, self.subgroups, self.moduni_db, self.T_K)

    # ── K-Values & Phase Boundaries ──────────────────────────────────────────

    def K_values(self, x: np.ndarray) -> np.ndarray:
        """Modified Raoult's law equilibrium ratios: K_i = gamma_i(x) * Psat_i(T) / P."""
        gamma = self._unifac.activity_coefficients(x)
        return np.array([
            gamma[i] * self.comps[i].vapor_pressure_Pa(self.T_C) / self.P_Pa
            for i in range(self.n)
        ])

    def bubble_point_pressure(self, x: np.ndarray) -> float:
        """Bubble-point pressure (Pa) at self.T_C for liquid composition x."""
        x = np.asarray(x, dtype=float)
        x = x / x.sum()
        gamma = self._unifac.activity_coefficients(x)
        return float(sum(
            x[i] * gamma[i] * self.comps[i].vapor_pressure_Pa(self.T_C)
            for i in range(self.n)
        ))

    def dew_point_pressure(self, y: np.ndarray, max_iter: int = 100, tol: float = 1e-10) -> float:
        """Iterative dew-point pressure (Pa) at self.T_C for vapour composition y."""
        y = np.asarray(y, dtype=float)
        y = y / y.sum()
        x = y.copy()
        for _ in range(max_iter):
            gamma = self._unifac.activity_coefficients(x)
            xn = np.array([
                y[i] / (gamma[i] * self.comps[i].vapor_pressure_Pa(self.T_C))
                for i in range(self.n)
            ])
            xn /= xn.sum()
            if np.max(np.abs(xn - x)) < tol:
                x = xn
                break
            x = xn
        gamma = self._unifac.activity_coefficients(x)
        return 1.0 / float(sum(
            y[i] / (gamma[i] * self.comps[i].vapor_pressure_Pa(self.T_C))
            for i in range(self.n)
        ))

    # ── Isothermal-Isobaric (PT) Flash ───────────────────────────────────────

    def flash(self, z: Union[List[float], np.ndarray], max_iter: int = 300, tol: float = 1e-10) -> dict:
        """Perform an isothermal-isobaric (PT) VLE flash at self.T_C, self.P_Pa."""
        z = np.array(z, dtype=float)
        z /= z.sum()

        P_bub = self.bubble_point_pressure(z)
        if self.P_Pa >= P_bub:
            return dict(
                phase="liquid", beta=0.0, x=z.copy(), y=np.zeros(self.n),
                K=self.K_values(z), gamma=self._unifac.activity_coefficients(z),
                P_bubble=P_bub, converged=True, iterations=0,
            )

        P_dew = self.dew_point_pressure(z)
        if self.P_Pa <= P_dew:
            return dict(
                phase="vapor", beta=1.0, x=np.zeros(self.n), y=z.copy(),
                K=self.K_values(z), gamma=self._unifac.activity_coefficients(z),
                P_bubble=P_bub, P_dew=P_dew, converged=True, iterations=0,
            )

        x = z.copy()
        beta, xn, yn = 0.5, z.copy(), z.copy()
        K = self.K_values(x)

        for it in range(max_iter):
            K = self.K_values(x)
            beta, xn, yn = rachford_rice(z, K)
            if np.max(np.abs(xn - x)) < tol:
                return dict(
                    phase="two-phase", beta=beta, x=xn, y=yn, K=K,
                    gamma=self._unifac.activity_coefficients(xn),
                    P_bubble=P_bub, P_dew=P_dew,
                    converged=True, iterations=it + 1,
                )
            x = xn

        return dict(
            phase="two-phase", beta=beta, x=xn, y=yn, K=K,
            gamma=self._unifac.activity_coefficients(xn),
            P_bubble=P_bub, P_dew=P_dew,
            converged=False, iterations=max_iter,
        )

    # ── Target Vapour-Fraction Pressure Finder ───────────────────────────────

    def find_pressure_for_target_beta(
        self,
        z: Union[List[float], np.ndarray],
        target_beta: float = 0.80,
        P_lo_mbar: float = 0.1,
        P_hi_mbar: Optional[float] = None,
        tol_beta: float = 1e-8,
        max_iter: int = 120,
    ) -> dict:
        """Bisect over pressure (mbar) to achieve exact `target_beta` at self.T_C."""
        z = np.array(z, dtype=float)
        z /= z.sum()

        if P_hi_mbar is None:
            P_hi_mbar = self.bubble_point_pressure(z) / 100.0 * 0.9999

        def beta_at(P_mb: float) -> Tuple[float, dict]:
            self.set_pressure(P_mb * 100.0)
            r = self.flash(z)
            return r["beta"], r

        b_hi, _ = beta_at(P_lo_mbar)
        b_lo, _ = beta_at(P_hi_mbar)

        if target_beta > b_hi:
            raise ValueError(
                f"Target beta={target_beta:.4f} exceeds maximum achievable "
                f"beta={b_hi:.6f} at P={P_lo_mbar} mbar."
            )
        if target_beta < b_lo:
            raise ValueError(
                f"Target beta={target_beta:.4f} is below beta={b_lo:.6f} "
                f"at P={P_hi_mbar:.4f} mbar."
            )

        lo, hi = P_lo_mbar, P_hi_mbar
        result = {}
        mid = 0.5 * (lo + hi)
        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            b_mid, result = beta_at(mid)
            if abs(b_mid - target_beta) < tol_beta:
                break
            if b_mid > target_beta:
                lo = mid
            else:
                hi = mid

        result["P_mbar"] = mid
        return result

    # ── Optional Post-Processing: Component Recovery Helper ──────────────────

    @staticmethod
    def compute_recovery(
        z: np.ndarray,
        flash_result: dict,
        feed_moles: np.ndarray,
    ) -> dict:
        """Optional helper to compute component split fractions and % recoveries."""
        beta = flash_result["beta"]
        x = flash_result["x"]
        y = flash_result["y"]
        n = len(z)

        F_total = float(feed_moles.sum())
        V_total = beta * F_total
        L_total = (1.0 - beta) * F_total

        alpha = np.zeros(n)
        moles_vap = np.zeros(n)
        moles_liq = np.zeros(n)
        R_vap = np.zeros(n)
        R_liq = np.zeros(n)

        for i in range(n):
            moles_vap[i] = V_total * y[i]
            moles_liq[i] = L_total * x[i]
            if z[i] > 1e-12:
                alpha[i] = beta * y[i] / z[i]
                R_vap[i] = alpha[i] * 100.0
                R_liq[i] = (1.0 - alpha[i]) * 100.0
            else:
                alpha[i] = R_vap[i] = R_liq[i] = float("nan")

        return dict(
            alpha=alpha,
            R_vap=R_vap,
            R_liq=R_liq,
            moles_vap=moles_vap,
            moles_liq=moles_liq,
            F_total=F_total,
            V_total=V_total,
            L_total=L_total,
        )

    def pure_component_table(self) -> List[dict]:
        """Return pure-component properties at self.T_C."""
        rows = []
        for k, c in zip(self.keys, self.comps):
            Psat = c.vapor_pressure_Pa(self.T_C)
            Tbp = c.saturation_temperature_C(101325.0)
            rows.append(dict(
                key=k, name=c.name, MW=c.MW,
                Tbp_C=Tbp, Psat_Pa=Psat, Psat_mbar=Psat / 100.0,
            ))
        return rows


# ==============================================================================
# STEP 2 VERIFICATION RUN
# ==============================================================================
if __name__ == "__main__":
    moduni_db = ModUNIFACDatabase("MODUNI.xlsx")

    requested = {
        "toluene": "Toluene",
        "water": "Water",
        "methanol": "Methanol",
        "DIPE": "Diisopropyl Ether",
        "EAc": "Ethyl acetate",
    }
    custom_antoine = {
        "DIPE": PureComponent(
            name="Diisopropyl Ether",
            MW=102.175,
            antoine_ABC=(6.86282, 1141.280, 215.570),
        )
    }
    comp_db = load_pure_components(
        requested_names=requested,
        perry_excel_path="PerryData.xlsx",
        custom_components=custom_antoine,
    )

    subgroups_map = {
        "toluene": {"ACH": 5, "ACCH3": 1},
        "water": {"H2O": 1},
        "methanol": {"CH3OH": 1},
        "DIPE": {"CH3": 4, "CH": 1, (13, "CHO"): 1},
        "EAc": {"CH3": 2, "CH2": 1, "CH3COO": 1},
    }

    keys = ["toluene", "water", "methanol", "DIPE", "EAc"]
    z_feed = [0.2, 0.2, 0.2, 0.2, 0.2]

    solver = FlashSolver(keys, comp_db, subgroups_map, moduni_db, T_C=60.0, P_Pa=70000.0)
    res = solver.flash(z_feed)

    print("Step 2 Check Passed! Flash at 60 °C, 70 kPa:")
    print(f"  Phase: {res['phase']} | Beta: {res['beta']:.4f} | P_bub: {res['P_bubble']/1000:.2f} kPa | P_dew: {res['P_dew']/1000:.2f} kPa | Iters: {res['iterations']}")
    for i, k in enumerate(keys):
        print(f"  {k:10s} | z={z_feed[i]:.3f} | x={res['x'][i]:.4f} | y={res['y'][i]:.4f} | K={res['K'][i]:.4f} | gamma={res['gamma'][i]:.4f}")

    vac = solver.find_pressure_for_target_beta(z_feed, target_beta=0.60)
    print(f"\nTarget beta=0.60 at 60 °C -> P = {vac['P_mbar']:.2f} mbar ({vac['P_mbar']/10:.2f} kPa)")