# Multicomponent Modified UNIFAC (Dortmund) VLE Flash Solver

An open-source Python thermodynamic engine and standalone Windows desktop GUI for multicomponent Vapor–Liquid Equilibrium (VLE) flash calculations, bubble/dew-point estimation, and target vapor-fraction pressure solving powered by the **Modified UNIFAC (Dortmund)** activity coefficient model.

> **Database Licensing Notice:** This repository distributes only open-source solver, verification, and GUI code under the MIT License. Proprietary or third-party thermodynamic parameter tables (`MODUNI.xlsx` and `PerryData.xlsx`) are **not** included. Users must supply their own Excel parameter files matching the schemas documented below, or run the solver using custom Antoine constants.

---

## Features

* **Modified UNIFAC (Dortmund) Engine (`VLE_flash_solver.py`):** Full implementation of the Weidlich & Gmehling (1987) / Gmehling et al. (1993) combinatorial ($3/4$-exponent volume fraction) and temperature-dependent residual group-contribution terms:
  $$\Psi_{nm}(T) = \exp\!\left(-\frac{a_{nm} + b_{nm}T + c_{nm}T^2}{T}\right)$$
* **Dual Calculation Modes:**
  * **PT Flash:** Given temperature $T$ and pressure $P$, determines phase state (`liquid`, `two-phase`, or `vapor`), vapor mole fraction $\beta$, equilibrium compositions $(x_i, y_i)$, $K$-values, and activity coefficients $\gamma_i$.
  * **Target Vapor Fraction ($\beta$) Pressure Finder:** Given $T$ and desired vapor split $\beta \in (0, 1)$, bisects over system pressure to solve for the exact operating pressure ($\text{kPa}$ and $\text{mbar}$).
* **Flexible Pure-Component Vapor Pressures:** Supports 5-parameter DIPPR 101 correlations via `PerryData.xlsx`, 3-parameter Antoine constants ( $\log_{10} P_{\text{mmHg}} = A - B/(C + T_{^\circ\text{C}})$ ), or a hybrid mix of both in the same mixture.
* **Standalone Windows GUI (`vle_gui.py`):** Zero-code desktop interface available as a precompiled `.exe` under [**Releases**](../../releases).
* **Automated Thermodynamic Verification (`verify_solver.py`):** Validates pure-component limits ($\gamma_i = 1$), numerical Gibbs–Duhem consistency, NRTL binary cross-checks, experimental azeotrope prediction (Ethanol–Water at $1\text{ atm}$), and flash mass-balance closure.

---

## Repository Structure

```text
├── VLE_flash_solver.py   # Core Modified UNIFAC + Rachford-Rice + FlashSolver engine
├── vle_gui.py            # Tkinter desktop GUI application
├── verify_solver.py      # Automated thermodynamic & azeotrope verification suite
├── moduni.ico            # Application icon bundled into the Windows executable
├── LICENSE               # MIT License
├── README.md             # Setup, Excel database schemas, and build instructions
└── .gitignore            # Excludes local .xlsx databases and PyInstaller build folders
```

---

## Quick Start: Non-Coders (Standalone Windows `.exe`)

1. Download **`VLE-Flash-Solver-v1.0.0.exe`** from the [**Releases**](../../releases) page.
2. Place your `MODUNI.xlsx` file (and optional `PerryData.xlsx`) in the same folder, or select them using the **Browse...** buttons inside the app.
3. Double-click **`VLE-Flash-Solver-v1.0.0.exe`** to launch the calculator—no Python installation required.

---

## Quick Start: Developers & Python Users

### 1. Installation
```bash
pip install numpy pandas scipy openpyxl
```

### 2. Launch the GUI from Source
```bash
python vle_gui.py
```

### 3. Run the Verification Suite
```bash
python verify_solver.py
```

### 4. Python Library Usage Example
```python
from VLE_flash_solver import PureComponent, load_pure_components, ModUNIFACDatabase, FlashSolver

# 1. Load Modified UNIFAC (Dortmund) subgroup & interaction tables
moduni_db = ModUNIFACDatabase("MODUNI.xlsx")

# 2. Load pure-component vapor pressure parameters (DIPPR 101 from PerryData + custom Antoine)
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

# 3. Define UNIFAC subgroup decompositions
subgroups = {
    "toluene": {"ACH": 5, "ACCH3": 1},
    "water": {"H2O": 1},
    "methanol": {"CH3OH": 1},
    "DIPE": {"CH3": 4, "CH": 1, (13, "CHO"): 1},  # Tuple resolves duplicate 'CHO' subgroup
    "EAc": {"CH3": 2, "CH2": 1, "CH3COO": 1},
}

# 4. Run PT Flash at 60 °C, 70 kPa
keys = ["toluene", "water", "methanol", "DIPE", "EAc"]
z_feed = [0.2, 0.2, 0.2, 0.2, 0.2]

solver = FlashSolver(keys, comp_db, subgroups, moduni_db, T_C=60.0, P_Pa=70000.0)
result = solver.flash(z_feed)

print(f"Phase: {result['phase']}, Vapor fraction (beta): {result['beta']:.4f}")
print("Liquid x:", result["x"])
print("Vapor  y:", result["y"])

# 5. Solve for operating pressure to achieve target beta = 0.60 at 60 °C
vac_result = solver.find_pressure_for_target_beta(z_feed, target_beta=0.60)
print(f"Required Pressure: {vac_result['P_mbar']:.2f} mbar")
```

---

## Required Excel Database Schemas

### 1. `MODUNI.xlsx` (Required)
Modified UNIFAC (Dortmund) parameters can be populated from published Dortmund parameter tables (Weidlich & Gmehling, *Ind. Eng. Chem. Res.*, 1987; Gmehling et al., 1993).

The workbook must contain **4 sheets**:

1. **`Groups` sheet** (Row 1 blank, Row 2 column headers, data starting at Row 3):
   | Column | Header | Type | Description |
   | :--- | :--- | :--- | :--- |
   | B | `Sr No` | `int` | Subgroup serial index |
   | C | `Subgrp` | `str` | Subgroup symbol (e.g., `CH3`, `ACH`, `H2O`, `CH3COO`) |
   | D | `Main grp` | `int` | Main group ID number ($1, 2, \dots, N$) |
   | E | `Main grp` | `str` | Main group name (e.g., `CH2`, `ACH`, `H2O`) |
   | F | `R` | `float` | Subgroup relative van der Waals volume parameter $R_k$ |
   | G | `Q` | `float` | Subgroup relative van der Waals surface area parameter $Q_k$ |

2. **`A`, `B`, and `C` sheets** (Main-group interaction parameter matrices):
   * **Row 2 (Columns C onward):** Integer Main Group IDs ($m = 1, 2, \dots, N$).
   * **Column B (Rows 3 onward):** Integer Main Group IDs ($n = 1, 2, \dots, N$).
   * **Matrix Cells $(n, m)$:** Interaction parameters $a_{nm}\text{ [K]}$, $b_{nm}\text{ [-]}$, and $c_{nm}\text{ [K}^{-1}\text{]}$ for $\Psi_{nm}(T)$.

#### Disambiguating Duplicate Subgroup Names (`CHO`)
In Dortmund parameter tables, `CHO` appears under two distinct main groups:
* **Main Group 10 (`CHO`):** Aldehyde `-CH=O` $\rightarrow$ specify as `(10, "CHO")`
* **Main Group 13 (`CH2O`):** Ether `>CH-O-` $\rightarrow$ specify as `(13, "CHO")`

All unique subgroup symbols are entered as plain strings (e.g., `"CH3"`, `"ACH"`, `"H2O"`).

---

### 2. `PerryData.xlsx` (Optional)
If `PerryData.xlsx` is omitted, set `perry_excel_path=None` in Python or click **Clear (Use Antoine Only)** in the GUI and enter 3-parameter Antoine constants.

When provided, `PerryData.xlsx` must contain:
1. **`Pvap` sheet** with columns `CmpdNo`, `Name`, `CASNo`, `C1`, `C2`, `C3`, `C4`, `C5`, `Tmin`, `Tmax K` for the DIPPR 101 equation ($T$ in $\text{K}$, $P^{\text{sat}}$ in $\text{Pa}$):
   $$\ln P^{\text{sat}}_{\text{Pa}} = C_1 + \frac{C_2}{T} + C_3 \ln T + C_4 T^{C_5}$$
2. **`LiqDens` sheet** with columns `CmpdNo` and `MolWt` ($\text{g/mol}$).

---

## Building the Standalone Windows Executable

To compile the GUI into a single Windows executable from source:

```bash
pip install pyinstaller pillow
pyinstaller --onefile --windowed --clean --icon="moduni.ico" --add-data "moduni.ico;." --name "VLE-Flash-Solver-v1.x.x" vle_gui.py
```

The standalone binary is written to `dist/VLE_Flash_Calculator.exe`.

---

## References

1. Weidlich, U., & Gmehling, J. (1987). A modified UNIFAC model. 1. Prediction of VLE, $h^E$, and $\gamma^\infty$. *Industrial & Engineering Chemistry Research*, 26(7), 1372–1381.
2. Gmehling, J., Li, J., & Schiller, M. (1993). A modified UNIFAC model. 2. Present parameter matrix and results for different thermodynamic properties. *Industrial & Engineering Chemistry Research*, 32(1), 178–193.
3. Michelsen, M. L. (1982). The isothermal flash problem. Part I. Stability. *Fluid Phase Equilibria*, 9(1), 1–19.

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
