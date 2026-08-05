"""
fit_magic_formula.py - fit a load-dependent Magic Formula to tire data.

FROM A NOTEBOOK (Jupyter, Spyder, IPython) -- use the function:

    from fit_magic_formula import fit_file
    fit_file("mydata.xlsx", mode="lateral", fix_c=1.0, linear_k=True)
    fit_file()                      # no argument -> file picker opens
    fit_file("sae_data.xlsx", flip_sign=True)   # if your Fy is negative

Do NOT run this file with %run inside a notebook. Jupyter puts its own arguments in sys.argv, the command-line parser reads them, fails, and raises SystemExit.

FROM A TERMINAL - run the file:

    python fit_magic_formula.py                      # file picker opens
    python fit_magic_formula.py mydata.xlsx --mode lateral
    python fit_magic_formula.py mydata.xlsx --mode longitudinal --shift

DATA FILE
Excel (.xlsx, .xls) or CSV. One header row, then three columns in this order. Column names are ignored; the order is what matters.

        slip        slip angle in DEGREES (lateral) or slip ratio (longitudinal)
        Fz          vertical load, N
        F           measured force, N

You need at least three different loads. The whole point is to fit the LOAD DEPENDENCE, and a single load cannot show it.

MODEL
    dfz = (Fz - Fz0) / Fz0
    mu  = pD1 + pD2*dfz
    D   = mu*Fz / sin(C*pi/2)   if C < 1, else mu*Fz
    C   = pC1
    E   = pE1 + pE2*dfz
    K   = pK1*Fz0*sin(2*arctan(Fz/(pK2*Fz0)))     or pK1*Fz with --linear-k
    B   = K / (C*D)

    F   = D*sin(C*arctan(B*x - E*(B*x - arctan(B*x)))) + SV

READING C
    C > 1   curve peaks then falls away. Normal on a hard surface.
    C ~ 1   rises to an asymptote, no falloff.
    C < 1   monotonic saturation, no peak. Typical of soft soil.
"""

import argparse
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

DEG = np.pi / 180.0
NAMES = ["pD1", "pD2", "pC1", "pE1", "pE2", "pK1", "pK2", "SV"]


def pick_file():
    """Open a file picker. Needs tkinter, which ships with standard Python."""
    import tkinter as tk
    from tkinter import filedialog
    root = tk.Tk()
    root.withdraw()
    path = filedialog.askopenfilename(
        title="Select tire data file",
        filetypes=[("Excel", "*.xlsx *.xls"), ("CSV", "*.csv"), ("All files", "*.*")])
    root.destroy()
    return path


def read_data(path):
    """First three columns are slip, Fz [N], F [N]. Names are ignored."""
    if path.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)
    df = df.iloc[:, :3].dropna()
    print(f"Read {path}")
    print(f"  columns used: {list(df.columns[:3])}")
    return df.iloc[:, 0].to_numpy(float), df.iloc[:, 1].to_numpy(float), df.iloc[:, 2].to_numpy(float)


def model(x, Fz, p, Fz0, linear_k):
    pD1, pD2, pC1, pE1, pE2, pK1, pK2, SV = p
    dfz = (Fz - Fz0) / Fz0

    mu = pD1 + pD2 * dfz
    C = pC1
    # mu is the largest force the curve actually reaches, divided by Fz.
    # For C < 1 the formula never gets to D, it tends to D*sin(C*pi/2), so
    # without this rescaling C and D trade off and the fit is meaningless.
    D = mu * Fz / (1.0 if C >= 1.0 else np.sin(C * np.pi / 2))
    E = pE1 + pE2 * dfz
    K = pK1 * Fz if linear_k else pK1 * Fz0 * np.sin(2.0 * np.arctan(Fz / (pK2 * Fz0)))
    B = K / (C * D)

    Bx = B * x
    return D * np.sin(C * np.arctan(Bx - E * (Bx - np.arctan(Bx)))) + SV * Fz / Fz0


def fit(x, Fz, F, Fz0, shift, fix_c, linear_k):
    #        pD1   pD2   pC1   pE1    pE2   pK1    pK2    SV
    lo = [0.01, -2.0, 0.45, -20.0, -5.0, 1e-3, 0.30, -500.0]
    hi = [5.00, 2.0, 2.20, 1.0, 5.0, 500.0, 12.0, 500.0]
    # Starting guess taken from the data. This matters: the residual is not convex, and starting pK1 at an arbitrary value drops the solver into a local minimum that fits just as well but returns a POSITIVE pD2, i.e. friction rising with load, which no tire does.
    near_origin = np.abs(x) < 0.25 * np.abs(x).max()
    slope = np.polyfit(x[near_origin], F[near_origin], 1)[0]   # N/rad at mean load

    if slope < 0:
        raise ValueError(
            "Force falls as slip rises, so your data is in the SAE sign "
            "convention (positive slip angle -> negative force). This model "
            "uses ISO. Pass flip_sign=True (or --flip-sign) to negate the "
            "force column, or negate it yourself in the spreadsheet.")

    # pK1 means different things in the two stiffness laws, so scale the starting guess to match whichever one is in use
    if linear_k:
        pK1_0 = slope / np.mean(Fz)
    else:
        pK1_0 = slope / (Fz0 * np.sin(2.0 * np.arctan(np.mean(Fz) / (1.9 * Fz0))))

    p0 = [np.max(np.abs(F)) / np.max(Fz), -0.05, 1.2, 0.0, 0.0, pK1_0, 1.9, 0.0]

    if not shift:
        lo[7] = hi[7] = p0[7] = 0.0
    if fix_c is not None:
        lo[2] = hi[2] = p0[2] = fix_c
    if linear_k:
        lo[6] = hi[6] = p0[6] = 1.0

    # bounds must be strictly ordered for the solver, so open pinned ones by eps
    lo = [v - 1e-9 for v in lo]
    hi = [v + 1e-9 for v in hi]

    # least_squares rejects an x0 outside the bounds, so bring it inside
    p0 = np.clip(p0, lo, hi)

    r = least_squares(lambda p: model(x, Fz, p, Fz0, linear_k) - F,
                      p0, bounds=(lo, hi), x_scale="jac", max_nfev=20000)

    pinned = {2 if fix_c is not None else -1, 6 if linear_k else -1, 7 if not shift else -1}
    at_bound = [NAMES[i] for i in range(8) if i not in pinned
                and (abs(r.x[i] - lo[i]) < 1e-6 or abs(r.x[i] - hi[i]) < 1e-6)]
    return r.x, at_bound


def fit_file(file=None, mode="lateral", fz0=None, shift=False,
             fix_c=None, linear_k=False, plot=None, flip_sign=False):
    """
    Fit a file and print the result. This is the notebook entry point:

        from fit_magic_formula import fit_file
        fit_file("mydata.xlsx", mode="lateral", fix_c=1.0, linear_k=True)

    Leave file as None to open a picker. Returns the coefficient array.
    """
    path = file if file else pick_file()
    slip, Fz, F = read_data(path)
    if flip_sign:
        F = -F                       # SAE data -> ISO convention
    x = slip * DEG if mode == "lateral" else slip
    loads = np.unique(Fz)
    Fz0 = fz0 if fz0 else float(np.median(loads))

    print(f"  {len(F)} points, {len(loads)} loads: {', '.join(f'{L:.0f}' for L in loads)} N")
    print(f"  mode {mode}, Fz0 = {Fz0:.1f} N")

    p, at_bound = fit(x, Fz, F, Fz0, shift, fix_c, linear_k)
    resid = model(x, Fz, p, Fz0, linear_k) - F
    rms = np.sqrt(np.mean(resid**2))
    r2 = 1.0 - np.sum(resid**2) / np.sum((F - F.mean())**2)

    print("\nFITTED COEFFICIENTS")
    for n, v in zip(NAMES, p):
        if n == "SV" and not shift:
            continue
        if n == "pK2" and linear_k:
            continue
        print(f"  {n:<4} = {v:+.5f}")
    print(f"\n  RMS {rms:.2f} N   max err {np.abs(resid).max():.1f} N   R2 {r2:.5f}")

    if at_bound:
        print(f"\n  AT BOUND: {', '.join(at_bound)}")
        print("  Your data does not determine these. Pin them with fix_c or")
        print("  simplify with linear_k, and say so when you report the fit.")

    s = "y" if mode == "lateral" else "x"
    print("\nPaste into on_road.py or off_road.py:")
    print(f"    Fz0={Fz0:.1f},")
    print(f"    pC{s}1={p[2]:.4f}, pD{s}1={p[0]:.4f}, pD{s}2={p[1]:.4f},")
    print(f"    pE{s}1={p[3]:.4f}, pE{s}2={p[4]:.4f}, pK{s}1={p[5]:.4f}, pK{s}2={p[6]:.4f},")
    if shift:
        print(f"    SV{s}={p[7]:.4f},")
    if linear_k:
        print("    linear_k=True,")

    if plot:
        import matplotlib.pyplot as plt
        conv = 1 / DEG if mode == "lateral" else 1.0
        fig, ax = plt.subplots(figsize=(7, 4.6))
        colors = plt.cm.viridis(np.linspace(0.05, 0.85, len(loads)))
        for c, L in zip(colors, loads):
            m = Fz == L
            ax.plot(x[m] * conv, F[m], "o", ms=4, color=c, label=f"{L:.0f} N")
            xs = np.linspace(0, x[m].max(), 300)
            ax.plot(xs * conv, model(xs, np.full_like(xs, L), p, Fz0, linear_k),
                    "-", color=c, lw=1.4)
        ax.set_xlabel("Slip angle [deg]" if mode == "lateral" else "Slip ratio [-]")
        ax.set_ylabel("Force [N]")
        ax.set_title(f"Magic Formula fit, C = {p[2]:.3f}")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, title="$F_z$")
        fig.tight_layout()
        fig.savefig(plot, dpi=200)
        print(f"\n  plot written to {plot}")

    return p


def main():
    ap = argparse.ArgumentParser(description="Fit a Magic Formula to tire data.")
    ap.add_argument("file", nargs="?", help="Excel or CSV file; omit to open a picker")
    ap.add_argument("--mode", choices=["lateral", "longitudinal"], default="lateral")
    ap.add_argument("--fz0", type=float, help="nominal load [N]; default is the median load")
    ap.add_argument("--shift", action="store_true", help="fit an offset SV (use for drawbar pull)")
    ap.add_argument("--fix-c", type=float, help="pin C instead of fitting it")
    ap.add_argument("--linear-k", action="store_true", help="use K = pK1*Fz")
    ap.add_argument("--flip-sign", action="store_true",
                    help="negate the force column (use for SAE-convention data)")
    ap.add_argument("--plot", help="write a fit plot to this path")
    a = ap.parse_args()

    fit_file(a.file, a.mode, a.fz0, a.shift, a.fix_c, a.linear_k,
             a.plot, a.flip_sign)


if __name__ == "__main__":
    main()