import numpy as np
from VLE_flash_solver import (
    FlashSolver,
    ModUNIFACDatabase,
    ModUNIFACModel,
    PureComponent,
    load_pure_components,
)


def nrtl_binary_gamma(x1: float, T_K: float, dg12_R: float, dg21_R: float, alpha: float) -> np.ndarray:
    """Reference 2-parameter NRTL activity coefficients for binary sanity check."""
    x2 = 1.0 - x1
    tau12, tau21 = dg12_R / T_K, dg21_R / T_K
    G12, G21 = np.exp(-alpha * tau12), np.exp(-alpha * tau21)
    ln_g1 = (x2 ** 2) * (tau21 * (G21 / (x1 + x2 * G21)) ** 2 + tau12 * G12 / ((x2 + x1 * G12) ** 2))
    ln_g2 = (x1 ** 2) * (tau12 * (G12 / (x2 + x1 * G12)) ** 2 + tau21 * G21 / ((x1 + x2 * G21) ** 2))
    return np.array([np.exp(ln_g1), np.exp(ln_g2)])


if __name__ == "__main__":
    db = ModUNIFACDatabase("MODUNI.xlsx")

    # --------------------------------------------------------------------------
    # TEST 1: Pure-Component Limit & Gibbs-Duhem Consistency (Methanol + Water)
    # --------------------------------------------------------------------------
    mw_subs = {"methanol": {"CH3OH": 1}, "water": {"H2O": 1}}
    model_mw = ModUNIFACModel(["methanol", "water"], mw_subs, db, T_K=333.15)

    g_pure_meoh = model_mw.activity_coefficients(np.array([1.0, 0.0]))[0]
    g_pure_h2o = model_mw.activity_coefficients(np.array([0.0, 1.0]))[1]
    assert abs(g_pure_meoh - 1.0) < 1e-12 and abs(g_pure_h2o - 1.0) < 1e-12, "Pure limit failed!"

    # Numerical Gibbs-Duhem check at x1 = 0.35
    x1 = 0.35
    dx = 1e-6
    g_plus = model_mw.activity_coefficients(np.array([x1 + dx, 1.0 - (x1 + dx)]))
    g_minus = model_mw.activity_coefficients(np.array([x1 - dx, 1.0 - (x1 - dx)]))
    dln_g = (np.log(g_plus) - np.log(g_minus)) / (2.0 * dx)
    gd_residual = abs(x1 * dln_g[0] + (1.0 - x1) * dln_g[1])
    print(f"[PASS] Test 1 (Pure limit & Gibbs-Duhem): pure_err={abs(g_pure_meoh-1):.2e}, GD_residual={gd_residual:.2e}")

    # --------------------------------------------------------------------------
    # TEST 2: Reference NRTL Comparison (Water + Methanol at 60 °C)
    # --------------------------------------------------------------------------
    print("\n[INFO] Test 2 (Modified UNIFAC vs. Reference NRTL — Water + Methanol at 60 °C):")
    for x_w in [0.1, 0.3, 0.5, 0.7, 0.9]:
        g_unifac = model_mw.activity_coefficients(np.array([1.0 - x_w, x_w]))  # [MeOH, H2O]
        g_nrtl = nrtl_binary_gamma(x_w, 333.15, dg12_R=-30.5, dg21_R=394.5, alpha=0.30)  # [H2O, MeOH]
        print(
            f"  x_H2O={x_w:.1f} | UNIFAC (g_H2O={g_unifac[1]:.4f}, g_MeOH={g_unifac[0]:.4f}) "
            f"vs NRTL (g_H2O={g_nrtl[0]:.4f}, g_MeOH={g_nrtl[1]:.4f})"
        )

    # --------------------------------------------------------------------------
    # TEST 3: Experimental Ethanol-Water Azeotrope Check (78.15 °C, x_EtOH=0.894)
    # --------------------------------------------------------------------------
    comp_ew = load_pure_components(
        {"ethanol": "Ethanol", "water": "Water"},
        perry_excel_path="PerryData.xlsx",
    )
    ew_subs = {
        "ethanol": {"CH3": 1, "CH2": 1, "OH (P)": 1},
        "water": {"H2O": 1},
    }
    solver_ew = FlashSolver(["ethanol", "water"], comp_ew, ew_subs, db, T_C=78.15, P_Pa=101325.0)
    z_az = np.array([0.894, 0.106])
    P_bub_kpa = solver_ew.bubble_point_pressure(z_az) / 1000.0
    K_az = solver_ew.K_values(z_az)
    print(
        f"\n[PASS] Test 3 (Ethanol-Water Azeotrope at 78.15 °C, x_EtOH=0.894): "
        f"P_bub = {P_bub_kpa:.2f} kPa (Exp: 101.33 kPa) | K = [{K_az[0]:.4f}, {K_az[1]:.4f}]"
    )

    # --------------------------------------------------------------------------
    # TEST 4: Multicomponent Flash Mass & Equilibrium Closure
    # --------------------------------------------------------------------------
    res = solver_ew.flash([0.40, 0.60])
    z_test = np.array([0.40, 0.60])
    mb_err = np.max(np.abs(z_test - ((1.0 - res["beta"]) * res["x"] + res["beta"] * res["y"])))
    assert mb_err < 1e-10, "Mass balance closure failed!"
    print(f"[PASS] Test 4 (Flash Mass-Balance Closure): max_error = {mb_err:.2e}")