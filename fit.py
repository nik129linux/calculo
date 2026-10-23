"""Least-squares line fit for a pendulum: T^2 = a L + b, so g = 4 pi^2 / a.

The fit is implemented from the hand-derived normal equations first; numpy and
scipy are used only afterwards, as an independent check.

Usage: python3 fit.py [datos.csv]   -> prints results, writes figs/ and resultados.md
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

G_REF = 9.78  # m/s^2, reference value for Colombia (approx.)


def load(path):
    """Rows: L_m, then one or more times for 10 oscillations. Returns L and T^2."""
    L, T2 = [], []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            times = [float(v) for k, v in row.items() if k.startswith("t10") and v]
            L.append(float(row["L_m"]))
            T2.append((sum(times) / len(times) / 10) ** 2)
    return np.array(L), np.array(T2)


def S(a, b, x, y):
    return float(np.sum((y - a * x - b) ** 2))


def gradient(a, b, x, y):
    # dS/da = -2 sum x_i (y_i - a x_i - b),  dS/db = -2 sum (y_i - a x_i - b)
    e = y - a * x - b
    return np.array([-2 * np.sum(x * e), -2 * np.sum(e)])


def fit_own(x, y):
    # grad S = 0 gives the normal equations:
    #   a Sxx + b Sx = Sxy
    #   a Sx  + b n  = Sy
    # solved by Cramer's rule, exactly as done by hand.
    n = x.size
    sx, sy = np.sum(x), np.sum(y)
    sxx, sxy = np.sum(x * x), np.sum(x * y)
    det = n * sxx - sx**2
    a = (n * sxy - sx * sy) / det
    b = (sxx * sy - sx * sxy) / det
    return float(a), float(b)


def hessian(x):
    # Constant: it does not depend on a, b or y, so S is a convex paraboloid.
    return np.array([[2 * np.sum(x * x), 2 * np.sum(x)], [2 * np.sum(x), 2 * x.size]])


def numeric_hessian(a, b, x, y, h=1e-3):
    # Central finite differences of S.
    saa = (S(a + h, b, x, y) - 2 * S(a, b, x, y) + S(a - h, b, x, y)) / h**2
    sbb = (S(a, b + h, x, y) - 2 * S(a, b, x, y) + S(a, b - h, x, y)) / h**2
    sab = (S(a + h, b + h, x, y) - S(a + h, b - h, x, y)
           - S(a - h, b + h, x, y) + S(a - h, b - h, x, y)) / (4 * h**2)
    return np.array([[saa, sab], [sab, sbb]])


def load_manual(path):
    # Hand-computed values, if present.
    p = Path(path)
    if not p.exists():
        return None
    with open(p) as f:
        return json.load(f)


def comparison_table(x, y, manual=None):
    # Markdown table with analytic, own and library values.
    from scipy import stats
    a, b = fit_own(x, y)
    s = S(a, b, x, y)
    pa, pb = np.polyfit(x, y, 1)
    lr = stats.linregress(x, y)
    s_poly = S(pa, pb, x, y)
    s_lr = S(lr.slope, lr.intercept, x, y)
    man = manual or {}
    rows = [
        ("a", man.get("a"), a, pa, lr.slope),
        ("b", man.get("b"), b, pb, lr.intercept),
        ("S_min", man.get("S_min"), s, s_poly, s_lr),
    ]
    out = ["| Parámetro | Analítico | Código propio | numpy.polyfit | scipy.linregress | Diferencia máxima |",
           "|---|---|---|---|---|---|"]
    for name, mval, own, npv, spv in rows:
        cell = f"{mval:.10f}" if mval is not None else "pendiente"
        vals = [v for v in (mval, own, npv, spv) if v is not None]
        diff = max(vals) - min(vals)
        out.append(f"| {name} | {cell} | {own:.10f} | {npv:.10f} | {spv:.10f} | {diff:.2e} |")
    return "\n".join(out)


def largest_residuals(a, b, x, y, k=3):
    # Indices of the k largest |residual|, largest first.
    e = np.abs(np.asarray(y) - a * np.asarray(x) - b)
    return np.argsort(e)[::-1][:k]


def quadratic_fit(x, y):
    # Quadratic least squares y = c2 x^2 + c1 x + c0.
    X = np.column_stack([np.asarray(x) ** 2, np.asarray(x), np.ones_like(np.asarray(x), dtype=float)])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    coef = np.asarray(coef, dtype=float)
    s = float(np.sum((np.asarray(y) - X @ coef) ** 2))
    cov = s / (x.size - 3) * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    return {"coef": coef, "S": s, "se": se}


def grid_ranges(a, b, x, y):
    # Symmetric box holding the S = 4 S_min ellipse, with margin.
    Hinv = np.linalg.inv(hessian(x))
    d = 3 * S(a, b, x, y)
    if abs(d) < 1e-12:
        d = 1e-6
    wa = 1.5 * math.sqrt(2 * d * Hinv[0, 0])
    wb = 1.5 * math.sqrt(2 * d * Hinv[1, 1])
    return (a - wa, a + wa, b - wb, b + wb)


def quality(a, b, x, y):
    e = y - a * x - b
    rmse = math.sqrt(np.mean(e**2))
    r2 = 1 - np.sum(e**2) / np.sum((y - y.mean()) ** 2)
    return rmse, float(r2)


def g_from_slope(a):
    return 4 * math.pi**2 / a


def sums(x, y):
    """The sums needed by hand, so the manual column of the table can be checked."""
    return {"n": x.size, "Sx": float(np.sum(x)), "Sy": float(np.sum(y)),
            "Sxx": float(np.sum(x * x)), "Sxy": float(np.sum(x * y))}


def stderr(a, b, x, y):
    # Residual variance s^2 = S_min / (n - 2), propagated to a and b
    n = x.size
    s2 = S(a, b, x, y) / (n - 2)
    sxx_c = np.sum((x - x.mean()) ** 2)
    return math.sqrt(s2 / sxx_c), math.sqrt(s2 * np.sum(x * x) / (n * sxx_c))


def length_offset(a, b):
    # Theory says b = 0. If T^2 = a (L + d), then b = a d: d is what the
    # measured lengths are missing (e.g. not measuring to the bob's centre).
    return b / a


def principal_directions(x):
    # Eigenvectors of the Hessian = axes of the elliptical level curves;
    # eigenvalues = curvature of S along each axis.
    return np.linalg.eigh(hessian(x))


def leave_one_out(x, y):
    """Option B: refit without each point; returns (i, a, b, g) per removed point."""
    out = []
    for i in range(x.size):
        m = np.arange(x.size) != i
        a, b = fit_own(x[m], y[m])
        out.append((i, a, b, g_from_slope(a)))
    return out


def plots(x, y, a, b, outdir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir.mkdir(exist_ok=True)
    e = y - a * x - b
    xs = np.linspace(0, x.max() * 1.05, 100)

    # 0. Scatter of the data only
    fig, ax = plt.subplots()
    ax.scatter(x, y, color="k")
    ax.set(xlabel="L (m)", ylabel="T² (s²)", title="Diagrama de dispersión")
    fig.savefig(outdir / "0_dispersion.png", dpi=150)

    # 1. Scatter + best line + worse lines
    fig, ax = plt.subplots()
    ax.scatter(x, y, color="k", zorder=3, label="datos")
    ax.plot(xs, a * xs + b, "C0", lw=2, label=f"óptima: S={S(a, b, x, y):.4f}")
    for da, db, c in [(0.4, -0.1, "C1"), (-0.4, 0.15, "C2")]:
        ax.plot(xs, (a + da) * xs + b + db, c, ls="--",
                label=f"a={a + da:.2f}, b={b + db:.2f}: S={S(a + da, b + db, x, y):.4f}")
    ax.set(xlabel="L (m)", ylabel="T² (s²)", title="Datos y rectas de ajuste")
    ax.legend(fontsize=8)
    fig.savefig(outdir / "1_dispersion_rectas.png", dpi=150)

    # Grid around the optimum for surface, traces and level curves
    a_lo, a_hi, b_lo, b_hi = grid_ranges(a, b, x, y)
    A, B = np.meshgrid(np.linspace(a_lo, a_hi, 200), np.linspace(b_lo, b_hi, 200))
    Z = np.sum((y[:, None, None] - A * x[:, None, None] - B) ** 2, axis=0)

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot_surface(A, B, Z, cmap="viridis", alpha=0.8)
    ax.scatter([a], [b], [S(a, b, x, y)], color="r", s=40)
    ax.set(xlabel="a", ylabel="b", zlabel="S(a,b)", title="Superficie de error")
    fig.savefig(outdir / "2_superficie.png", dpi=150)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    aa = np.linspace(a_lo, a_hi, 200)
    bb = np.linspace(b_lo, b_hi, 200)
    wa, wb = a_hi - a, b_hi - b  # half-widths of the window
    for db in (0, 0.5 * wb, -0.5 * wb):
        ax1.plot(aa, [S(t, b + db, x, y) for t in aa], label=f"b = {b + db:.3f}")
    for da in (0, 0.5 * wa, -0.5 * wa):
        ax2.plot(bb, [S(a + da, t, x, y) for t in bb], label=f"a = {a + da:.3f}")
    ax1.axvline(a, color="r", ls=":")
    ax2.axvline(b, color="r", ls=":")
    ax1.set(xlabel="a", ylabel="S", title="Trazas con b fijo")
    ax2.set(xlabel="b", ylabel="S", title="Trazas con a fijo")
    ax1.legend(fontsize=8)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "3_trazas.png", dpi=150)

    fig, ax = plt.subplots()
    cs = ax.contour(A, B, Z, levels=12, cmap="viridis")
    ax.clabel(cs, fontsize=7)
    ax.plot(a, b, "r*", ms=12, label=f"(a*, b*) = ({a:.3f}, {b:.3f})")
    lams, vecs = principal_directions(x)
    for lam, v, c in zip(lams, vecs.T, ("C3", "C1")):
        r = 0.8 * min(wa / abs(v[0]) if v[0] else np.inf, wb / abs(v[1]) if v[1] else np.inf)
        ax.plot([a - r * v[0], a + r * v[0]], [b - r * v[1], b + r * v[1]], c,
                label=f"autovector λ = {lam:.2f}")
    ax.set(xlabel="a", ylabel="b", title="Curvas de nivel de S(a,b)")
    ax.legend(fontsize=8)
    fig.savefig(outdir / "4_curvas_nivel.png", dpi=150)

    # Cross-sections through the minimum along the Hessian's eigenvectors
    fig, ax = plt.subplots()
    t = np.linspace(-math.hypot(wa, wb), math.hypot(wa, wb), 200)
    for lam, v, c in zip(lams, vecs.T, ("C3", "C1")):
        ax.plot(t, [S(a + s * v[0], b + s * v[1], x, y) for s in t], c,
                label=f"dirección ({v[0]:.2f}, {v[1]:.2f}), curvatura λ = {lam:.2f}")
    ax.set(xlabel="t (distancia desde (a*, b*))", ylabel="S",
           title="Secciones transversales en las direcciones principales")
    ax.legend(fontsize=8)
    fig.savefig(outdir / "6_secciones_principales.png", dpi=150)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.axhline(0, color="gray")
    ax1.scatter(x, e, color="k")
    ax1.set(xlabel="L (m)", ylabel="residuo e (s²)", title="Residuos vs L")
    ax2.hist(e, bins="auto", edgecolor="k")
    ax2.set(xlabel="residuo e (s²)", ylabel="frecuencia", title="Histograma de residuos")
    fig.tight_layout()
    fig.savefig(outdir / "5_residuos.png", dpi=150)

    plt.close("all")


def main(path):
    x, y = load(path)
    a, b = fit_own(x, y)
    s = S(a, b, x, y)
    H = hessian(x)
    Hn = numeric_hessian(a, b, x, y)
    hdiff = float(np.max(np.abs(Hn - H)))
    eig = np.linalg.eigvalsh(H)
    rmse, r2 = quality(a, b, x, y)
    g = g_from_slope(a)
    q = quadratic_fit(x, y)
    c2, c1, c0 = q["coef"]
    s_quad = q["S"]
    se_q = q["se"]
    manual = load_manual(Path(path).parent / "manual.json")
    table = comparison_table(x, y, manual)

    sm = sums(x, y)
    se_a, se_b = stderr(a, b, x, y)
    e = y - a * x - b
    idx = largest_residuals(a, b, x, y, k=3)
    lines = [
        f"n = {x.size} puntos",
        "",
        "Sumas para verificar el cálculo a mano:",
        "| i | x = L (m) | y = T² (s²) | x² | xy |",
        "|---|---|---|---|---|",
    ]
    for i, (xi, yi) in enumerate(zip(x, y)):
        lines.append(f"| {i} | {xi:.4f} | {yi:.5f} | {xi * xi:.5f} | {xi * yi:.5f} |")
    lines += [
        f"| Σ | {sm['Sx']:.4f} | {sm['Sy']:.5f} | {sm['Sxx']:.5f} | {sm['Sxy']:.5f} |",
        "",
    ]
    lines += table.splitlines()
    lines += [
        "",
        f"gradiente en (a*, b*) = {gradient(a, b, x, y)}",
        f"Hessiano = {H.tolist()}",
        f"Hessiano numérico = {Hn.tolist()}, diferencia máxima vs analítico = {hdiff:.2e}",
        f"S_aa = {H[0, 0]:.4f} > 0, det H = {np.linalg.det(H):.4f} > 0, autovalores = {eig}",
        f"RMSE = {rmse:.5f} s², R² = {r2:.6f}",
        f"g = 4π²/a = {g:.4f} m/s² (ref {G_REF}, error {abs(g - G_REF) / G_REF:.2%})",
        f"error estándar: a = {a:.4f} ± {se_a:.4f}, b = {b:.4f} ± {se_b:.4f}",
        f"g = {g:.3f} ± {g * se_a / a:.3f} m/s²",
        f"b / error(b) = {b / se_b:.2f} (|valor| > 2 → b distinto de 0, la teoría dice 0)",
        f"corrimiento de L implícito en b: d = b/a = {length_offset(a, b) * 100:.2f} cm",
        "",
        "Residuos mayores:",
    ]
    for i in idx:
        i = int(i)
        lines.append(f"  punto {i} (L={x[i]:.3f}): e={e[i]:.5f}")
    lines += [
        "",
        "Opción B (quitar un punto):",
    ]
    for i, ai, bi, gi in leave_one_out(x, y):
        lines.append(f"  sin punto {i} (L={x[i]:.3f}): a={ai:.4f}, b={bi:.4f}, g={gi:.4f}")
    lines += [
        "",
        f"Opción C: cuadrático y = {c2:.4f}x² + {c1:.4f}x + {c0:.4f}, S = {s_quad:.6f} vs lineal S = {s:.6f}",
        f"Opción C: c2 = {c2:.4f} ± {se_q[0]:.4f}, c2/se = {c2 / se_q[0]:.2f} (|valor| > 2 → la curvatura es significativa)",
    ]
    text = "\n".join(lines)
    print(text)
    out = Path(path).parent
    (out / "resultados.md").write_text(text + "\n")
    plots(x, y, a, b, out / "figs")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent / "datos.csv")
