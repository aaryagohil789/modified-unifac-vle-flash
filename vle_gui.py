import ast
import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import numpy as np

from VLE_flash_solver import (
    FlashSolver,
    ModUNIFACDatabase,
    PureComponent,
    load_pure_components,
)


def _resource_path(relative_path: str) -> str:
    """Resolve bundled asset paths for both source execution and PyInstaller --onefile."""
    base_path = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)


class VLEFlashGUI(tk.Tk):
    """
    Standalone Tkinter Desktop GUI for Modified UNIFAC (Dortmund) VLE Flashing.
    Designed to import `VLE_flash_solver.py` cleanly and compile into a single
    executable via PyInstaller for GitHub Releases.
    """

    def __init__(self):
        super().__init__()
        self.title("Modified UNIFAC (Dortmund) Multicomponent VLE Flash Calculator")
        self.geometry("1040x720")
        self.minsize(920, 640)

        try:
            self.iconbitmap(_resource_path("moduni.ico"))
        except Exception:
            pass

        self.moduni_db = None
        self._build_ui()
        self._load_default_demo_rows()

    # ── UI Construction ──────────────────────────────────────────────────────

    def _build_ui(self):
        pad = {"padx": 10, "pady": 6}

        # 1. Database File Selection Frame
        db_frame = ttk.LabelFrame(self, text="1. Thermodynamic Databases (User-Supplied Excel Files)")
        db_frame.pack(fill="x", **pad)

        ttk.Label(db_frame, text="MODUNI.xlsx (Required):").grid(row=0, column=0, sticky="w", padx=6, pady=4)
        self.moduni_path_var = tk.StringVar(value="MODUNI.xlsx")
        ttk.Entry(db_frame, textvariable=self.moduni_path_var, width=65).grid(row=0, column=1, padx=4, pady=4)
        ttk.Button(db_frame, text="Browse...", command=self._browse_moduni).grid(row=0, column=2, padx=4, pady=4)

        ttk.Label(db_frame, text="PerryData.xlsx (Optional):").grid(row=1, column=0, sticky="w", padx=6, pady=4)
        self.perry_path_var = tk.StringVar(value="PerryData.xlsx")
        ttk.Entry(db_frame, textvariable=self.perry_path_var, width=65).grid(row=1, column=1, padx=4, pady=4)
        ttk.Button(db_frame, text="Browse...", command=self._browse_perry).grid(row=1, column=2, padx=4, pady=4)
        ttk.Button(db_frame, text="Clear (Use Antoine Only)", command=lambda: self.perry_path_var.set("")).grid(
            row=1, column=3, padx=4, pady=4
        )

        # 2. Component Table Frame
        comp_frame = ttk.LabelFrame(
            self,
            text="2. Mixture Components, Vapor Pressure Source, & UNIFAC Subgroups",
        )
        comp_frame.pack(fill="both", expand=True, **pad)

        # Toolbar for adding/removing rows
        btn_bar = ttk.Frame(comp_frame)
        btn_bar.pack(fill="x", padx=6, pady=4)
        ttk.Button(btn_bar, text="+ Add Component", command=self._add_component_row).pack(side="left", padx=3)
        ttk.Button(btn_bar, text="- Remove Selected", command=self._remove_selected_row).pack(side="left", padx=3)
        ttk.Button(btn_bar, text="Load 5-Component Demo", command=self._load_default_demo_rows).pack(side="left", padx=12)

        hint = (
            "Hint: For 'CHO', specify tuple with main-group number: (10, 'CHO') = aldehyde, (13, 'CHO') = ether."
        )
        ttk.Label(btn_bar, text=hint, foreground="#444444").pack(side="right", padx=6)

        # Scrollable Component Editor Canvas
        canvas_container = ttk.Frame(comp_frame)
        canvas_container.pack(fill="both", expand=True, padx=6, pady=4)

        self.comp_canvas = tk.Canvas(canvas_container, height=190, highlightthickness=0)
        scrollbar = ttk.Scrollbar(canvas_container, orient="vertical", command=self.comp_canvas.yview)
        self.rows_frame = ttk.Frame(self.comp_canvas)

        self.rows_frame.bind(
            "<Configure>",
            lambda e: self.comp_canvas.configure(scrollregion=self.comp_canvas.bbox("all")),
        )
        self.comp_canvas.create_window((0, 0), window=self.rows_frame, anchor="nw")
        self.comp_canvas.configure(yscrollcommand=scrollbar.set)

        self.comp_canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Header Row
        headers = [
            ("Sel", 4),
            ("Name / Key", 16),
            ("Feed z_i", 10),
            ("Pvap Mode", 11),
            ("Antoine (A, B, C) [mmHg, °C]", 26),
            ("UNIFAC Subgroups Dict", 36),
        ]
        for col_idx, (text, w) in enumerate(headers):
            ttk.Label(self.rows_frame, text=text, width=w, font=("Segoe UI", 9, "bold")).grid(
                row=0, column=col_idx, padx=3, pady=3, sticky="w"
            )

        self.comp_rows = []

        # 3. Operating Conditions & Solver Mode Frame
        op_frame = ttk.LabelFrame(self, text="3. Operating Conditions & Calculation Mode")
        op_frame.pack(fill="x", **pad)

        self.calc_mode_var = tk.StringVar(value="pt_flash")
        ttk.Radiobutton(
            op_frame,
            text="PT Flash (Given T & P)",
            variable=self.calc_mode_var,
            value="pt_flash",
            command=self._toggle_mode_inputs,
        ).grid(row=0, column=0, padx=8, pady=5, sticky="w")
        ttk.Radiobutton(
            op_frame,
            text="Find Pressure for Target Vapor Fraction (β)",
            variable=self.calc_mode_var,
            value="target_beta",
            command=self._toggle_mode_inputs,
        ).grid(row=0, column=1, padx=8, pady=5, sticky="w")

        ttk.Label(op_frame, text="Temperature (°C):").grid(row=1, column=0, padx=8, pady=5, sticky="e")
        self.temp_var = tk.StringVar(value="60.0")
        ttk.Entry(op_frame, textvariable=self.temp_var, width=14).grid(row=1, column=1, padx=4, pady=5, sticky="w")

        ttk.Label(op_frame, text="Pressure (kPa):").grid(row=1, column=2, padx=8, pady=5, sticky="e")
        self.press_kpa_var = tk.StringVar(value="70.0")
        self.press_entry = ttk.Entry(op_frame, textvariable=self.press_kpa_var, width=14)
        self.press_entry.grid(row=1, column=3, padx=4, pady=5, sticky="w")

        ttk.Label(op_frame, text="Target β (0–1):").grid(row=1, column=4, padx=8, pady=5, sticky="e")
        self.target_beta_var = tk.StringVar(value="0.60")
        self.beta_entry = ttk.Entry(op_frame, textvariable=self.target_beta_var, width=14, state="disabled")
        self.beta_entry.grid(row=1, column=5, padx=4, pady=5, sticky="w")

        run_btn = ttk.Button(op_frame, text="▶ Run Flash Calculation", command=self._run_calculation)
        run_btn.grid(row=0, column=4, columnspan=2, padx=10, pady=5, sticky="ew")

        # 4. Results Output Frame
        res_frame = ttk.LabelFrame(self, text="4. Equilibrium Results")
        res_frame.pack(fill="both", expand=True, **pad)

        self.summary_var = tk.StringVar(value="Ready. Load databases and click 'Run Flash Calculation'.")
        ttk.Label(res_frame, textvariable=self.summary_var, font=("Segoe UI", 10, "bold")).pack(
            anchor="w", padx=8, pady=4
        )

        cols = ("comp", "z", "x", "y", "K", "gamma", "psat")
        self.tree = ttk.Treeview(res_frame, columns=cols, show="headings", height=7)
        headings = {
            "comp": "Component",
            "z": "Feed z_i",
            "x": "Liquid x_i",
            "y": "Vapor y_i",
            "K": "K-Value",
            "gamma": "Activity γ_i",
            "psat": "Psat (kPa)",
        }
        for c, label in headings.items():
            self.tree.heading(c, text=label)
            self.tree.column(c, width=125, anchor="center")
        self.tree.pack(fill="both", expand=True, padx=6, pady=4)

    # ── Mode Input Toggling ──────────────────────────────────────────────────

    def _toggle_mode_inputs(self):
        if self.calc_mode_var.get() == "pt_flash":
            self.press_entry.config(state="normal")
            self.beta_entry.config(state="disabled")
        else:
            self.press_entry.config(state="disabled")
            self.beta_entry.config(state="normal")

    # ── Row Management ───────────────────────────────────────────────────────

    def _add_component_row(self, name="", z="0.2", mode="PerryData", antoine="", subgroups=""):
        row_idx = len(self.comp_rows) + 1
        sel_var = tk.BooleanVar(value=False)
        name_var = tk.StringVar(value=name)
        z_var = tk.StringVar(value=z)
        mode_var = tk.StringVar(value=mode)
        ant_var = tk.StringVar(value=antoine)
        sg_var = tk.StringVar(value=subgroups)

        chk = ttk.Checkbutton(self.rows_frame, variable=sel_var)
        e_name = ttk.Entry(self.rows_frame, textvariable=name_var, width=18)
        e_z = ttk.Entry(self.rows_frame, textvariable=z_var, width=10)
        cb_mode = ttk.Combobox(
            self.rows_frame, textvariable=mode_var, values=["PerryData", "Antoine"], width=10, state="readonly"
        )
        e_ant = ttk.Entry(self.rows_frame, textvariable=ant_var, width=28)
        e_sg = ttk.Entry(self.rows_frame, textvariable=sg_var, width=40)

        widgets = [chk, e_name, e_z, cb_mode, e_ant, e_sg]
        for c_i, w in enumerate(widgets):
            w.grid(row=row_idx, column=c_i, padx=3, pady=2, sticky="w")

        self.comp_rows.append({
            "sel": sel_var,
            "name": name_var,
            "z": z_var,
            "mode": mode_var,
            "antoine": ant_var,
            "subgroups": sg_var,
            "widgets": widgets,
        })

    def _remove_selected_row(self):
        remaining = []
        for r in self.comp_rows:
            if r["sel"].get():
                for w in r["widgets"]:
                    w.destroy()
            else:
                remaining.append(r)
        self.comp_rows = remaining
        for idx, r in enumerate(self.comp_rows, start=1):
            for c_i, w in enumerate(r["widgets"]):
                w.grid(row=idx, column=c_i, padx=3, pady=2, sticky="w")

    def _load_default_demo_rows(self):
        for r in self.comp_rows:
            for w in r["widgets"]:
                w.destroy()
        self.comp_rows.clear()

        demo = [
            ("Toluene", "0.20", "PerryData", "", "{'ACH': 5, 'ACCH3': 1}"),
            ("Water", "0.20", "PerryData", "", "{'H2O': 1}"),
            ("Methanol", "0.20", "PerryData", "", "{'CH3OH': 1}"),
            ("Diisopropyl Ether", "0.20", "Antoine", "6.86282, 1141.280, 215.570", "{'CH3': 4, 'CH': 1, (13, 'CHO'): 1}"),
            ("Ethyl acetate", "0.20", "PerryData", "", "{'CH3': 2, 'CH2': 1, 'CH3COO': 1}"),
        ]
        for item in demo:
            self._add_component_row(*item)

    # ── File Browsers ────────────────────────────────────────────────────────

    def _browse_moduni(self):
        path = filedialog.askopenfilename(title="Select MODUNI.xlsx", filetypes=[("Excel files", "*.xlsx *.xls")])
        if path:
            self.moduni_path_var.set(path)
            self.moduni_db = None

    def _browse_perry(self):
        path = filedialog.askopenfilename(title="Select PerryData.xlsx", filetypes=[("Excel files", "*.xlsx *.xls")])
        if path:
            self.perry_path_var.set(path)

    # ── Execution Logic ──────────────────────────────────────────────────────

    def _run_calculation(self):
        try:
            moduni_path = self.moduni_path_var.get().strip()
            if not moduni_path:
                raise ValueError("Please select your MODUNI.xlsx database file.")

            if self.moduni_db is None or self.moduni_db.excel_path != moduni_path:
                self.moduni_db = ModUNIFACDatabase(moduni_path)

            perry_path = self.perry_path_var.get().strip() or None

            if not self.comp_rows:
                raise ValueError("Please add at least two components.")

            keys = []
            requested_names = {}
            custom_components = {}
            subgroups_map = {}
            z_list = []

            for r in self.comp_rows:
                name = r["name"].get().strip()
                if not name:
                    raise ValueError("All component rows must have a Name / Key.")
                keys.append(name)
                requested_names[name] = name
                z_list.append(float(r["z"].get().strip()))

                sg_raw = r["subgroups"].get().strip()
                sg_dict = ast.literal_eval(sg_raw)
                if not isinstance(sg_dict, dict) or not sg_dict:
                    raise ValueError(f"Invalid UNIFAC subgroup dictionary for '{name}'.")
                subgroups_map[name] = sg_dict

                if r["mode"].get() == "Antoine":
                    parts = [float(v.strip()) for v in r["antoine"].get().split(",")]
                    if len(parts) != 3:
                        raise ValueError(f"Enter exactly 3 comma-separated Antoine constants (A, B, C) for '{name}'.")
                    custom_components[name] = PureComponent(
                        name=name,
                        MW=0.0,
                        antoine_ABC=(parts[0], parts[1], parts[2]),
                    )

            comp_db = load_pure_components(
                requested_names=requested_names,
                perry_excel_path=perry_path,
                custom_components=custom_components,
            )

            T_C = float(self.temp_var.get().strip())
            z_arr = np.array(z_list, dtype=float)
            z_arr /= z_arr.sum()

            if self.calc_mode_var.get() == "pt_flash":
                P_kPa = float(self.press_kpa_var.get().strip())
                solver = FlashSolver(keys, comp_db, subgroups_map, self.moduni_db, T_C=T_C, P_Pa=P_kPa * 1000.0)
                res = solver.flash(z_arr)
                p_used_kpa = P_kPa
            else:
                target_b = float(self.target_beta_var.get().strip())
                solver = FlashSolver(keys, comp_db, subgroups_map, self.moduni_db, T_C=T_C, P_Pa=101325.0)
                res = solver.find_pressure_for_target_beta(z_arr, target_beta=target_b)
                p_used_kpa = res["P_mbar"] / 10.0

            # Populate Results Table
            for row in self.tree.get_children():
                self.tree.delete(row)

            pbub_kpa = res.get("P_bubble", float("nan")) / 1000.0
            pdew_kpa = res.get("P_dew", float("nan")) / 1000.0
            self.summary_var.set(
                f"Phase: {res['phase'].upper()}  |  Vapor Fraction (β): {res['beta']:.4f}  |  "
                f"P: {p_used_kpa:.3f} kPa ({p_used_kpa*10:.2f} mbar)  |  "
                f"P_bub: {pbub_kpa:.2f} kPa  |  P_dew: {pdew_kpa:.2f} kPa"
            )

            for i, k in enumerate(keys):
                psat_kpa = comp_db[k].vapor_pressure_Pa(T_C) / 1000.0
                self.tree.insert(
                    "",
                    "end",
                    values=(
                        k,
                        f"{z_arr[i]:.4f}",
                        f"{res['x'][i]:.4f}",
                        f"{res['y'][i]:.4f}",
                        f"{res['K'][i]:.4f}",
                        f"{res['gamma'][i]:.4f}",
                        f"{psat_kpa:.3f}",
                    ),
                )

        except Exception as exc:
            messagebox.showerror("Calculation Error", str(exc))


if __name__ == "__main__":
    app = VLEFlashGUI()
    app.mainloop()