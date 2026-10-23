import math

import numpy as np
import pytest
from scipy import stats

import fit

G = 9.78


def synthetic(noise=0.0, seed=0):
    # Ideal pendulum: T^2 = (4 pi^2 / g) L, plus optional noise on T^2
    rng = np.random.default_rng(seed)
    L = np.linspace(0.2, 1.2, 11)
    T2 = 4 * math.pi**2 / G * L + rng.normal(0, noise, L.size)
    return L, T2


def test_noise_free_recovers_g():
    x, y = synthetic()
    a, b = fit.fit_own(x, y)
    assert fit.g_from_slope(a) == pytest.approx(G, rel=1e-12)
    assert b == pytest.approx(0, abs=1e-12)


def test_own_matches_libraries():
    x, y = synthetic(noise=0.05)
    a, b = fit.fit_own(x, y)
    pa, pb = np.polyfit(x, y, 1)
    lr = stats.linregress(x, y)
    assert abs(a - pa) < 1e-9 and abs(b - pb) < 1e-9
    assert abs(a - lr.slope) < 1e-9 and abs(b - lr.intercept) < 1e-9


def test_gradient_zero_at_optimum():
    x, y = synthetic(noise=0.05)
    a, b = fit.fit_own(x, y)
    assert np.allclose(fit.gradient(a, b, x, y), 0, atol=1e-9)


def test_optimum_beats_neighbours():
    x, y = synthetic(noise=0.05)
    a, b = fit.fit_own(x, y)
    s = fit.S(a, b, x, y)
    for da, db in [(0.01, 0), (-0.01, 0), (0, 0.01), (0, -0.01), (0.01, -0.01)]:
        assert fit.S(a + da, b + db, x, y) > s


def test_hessian_matches_finite_differences_and_is_positive_definite():
    x, y = synthetic(noise=0.05)
    H = fit.hessian(x)
    a, b = fit.fit_own(x, y)
    h = 1e-4
    num_aa = (fit.S(a + h, b, x, y) - 2 * fit.S(a, b, x, y) + fit.S(a - h, b, x, y)) / h**2
    num_bb = (fit.S(a, b + h, x, y) - 2 * fit.S(a, b, x, y) + fit.S(a, b - h, x, y)) / h**2
    assert num_aa == pytest.approx(H[0, 0], rel=1e-5)
    assert num_bb == pytest.approx(H[1, 1], rel=1e-5)
    assert H[0, 0] > 0 and np.linalg.det(H) > 0
    assert np.all(np.linalg.eigvalsh(H) > 0)


def test_quality_metrics():
    x, y = synthetic()
    a, b = fit.fit_own(x, y)
    rmse, r2 = fit.quality(a, b, x, y)
    assert rmse == pytest.approx(0, abs=1e-12)
    assert r2 == pytest.approx(1, abs=1e-12)


def test_load_averages_trials(tmp_path):
    p = tmp_path / "d.csv"
    p.write_text("L_m,t10_1_s,t10_2_s,t10_3_s\n0.5,14.0,14.2,14.1\n1.0,20.0,20.0,20.0\n")
    x, y = fit.load(p)
    assert x.tolist() == [0.5, 1.0]
    assert y[0] == pytest.approx(1.41**2)
    assert y[1] == pytest.approx(4.0)


def test_sums_for_hand_check():
    x, y = synthetic(noise=0.05)
    s = fit.sums(x, y)
    assert s["n"] == x.size
    assert s["Sxy"] == pytest.approx(np.sum(x * y))
    assert s["Sxx"] == pytest.approx(np.sum(x * x))


def test_intercept_stderr_matches_scipy():
    x, y = synthetic(noise=0.05)
    a, b = fit.fit_own(x, y)
    lr = stats.linregress(x, y)
    se_a, se_b = fit.stderr(a, b, x, y)
    assert se_a == pytest.approx(lr.stderr, rel=1e-9)
    assert se_b == pytest.approx(lr.intercept_stderr, rel=1e-9)


def test_length_offset_from_intercept():
    # L measured 2 cm short: true T^2 = k (L + 0.02) -> b = k * 0.02, offset = b/a = +0.02
    x, _ = synthetic()
    k = 4 * math.pi**2 / G
    a, b = fit.fit_own(x, k * (x + 0.02))
    assert fit.length_offset(a, b) == pytest.approx(0.02)


def test_curvature_along_eigenvectors_equals_eigenvalues():
    x, y = synthetic(noise=0.05)
    a, b = fit.fit_own(x, y)
    lams, vecs = fit.principal_directions(x)
    h = 1e-3
    for lam, v in zip(lams, vecs.T):
        f = lambda t: fit.S(a + t * v[0], b + t * v[1], x, y)
        assert (f(h) - 2 * f(0) + f(-h)) / h**2 == pytest.approx(lam, rel=1e-5)


# --- Review fixes (sections 5, 7.6, 8.7, 9, 10, 13C of the assignment) ---

def test_numeric_hessian_matches_analytic():
    x, y = synthetic(noise=0.05)
    a, b = fit.fit_own(x, y)
    Hn = fit.numeric_hessian(a, b, x, y)
    assert Hn.shape == (2, 2)
    assert np.allclose(Hn, fit.hessian(x), rtol=1e-4)
    assert Hn[0, 1] == pytest.approx(Hn[1, 0], rel=1e-6)


def test_comparison_table_has_manual_column():
    x, y = synthetic(noise=0.05)
    t = fit.comparison_table(x, y)
    assert "Analítico" in t.splitlines()[0]
    assert t.count("pendiente") == 3  # a, b, S_min without hand values
    a, b = fit.fit_own(x, y)
    t2 = fit.comparison_table(x, y, manual={"a": round(a, 2), "b": round(b, 2), "S_min": 0.123})
    assert "pendiente" not in t2
    row_a = next(l for l in t2.splitlines() if l.startswith("| a "))
    # max difference must now include the hand-rounded value, so it is >= the rounding gap
    assert float(row_a.rstrip(" |").split("|")[-1]) >= abs(a - round(a, 2)) * 0.999


def test_load_manual(tmp_path):
    assert fit.load_manual(tmp_path / "missing.json") is None
    p = tmp_path / "manual.json"
    p.write_text('{"a": 4.01, "b": 0.02, "S_min": 0.011}')
    assert fit.load_manual(p) == {"a": 4.01, "b": 0.02, "S_min": 0.011}


def test_largest_residuals():
    x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    y = 2 * x + 1
    y[3] += 5.0
    y[1] -= 1.0
    a, b = fit.fit_own(x, y)
    idx = fit.largest_residuals(a, b, x, y, k=2)
    assert list(idx)[0] == 3 and len(idx) == 2


def test_quadratic_fit_stderr():
    x, y = synthetic(noise=0.05, seed=3)
    q = fit.quadratic_fit(x, y)
    c, s, se = q["coef"], q["S"], q["se"]
    X = np.column_stack([x**2, x, np.ones_like(x)])
    ref, *_ = np.linalg.lstsq(X, y, rcond=None)
    assert np.allclose(c, ref)
    assert s == pytest.approx(float(np.sum((y - X @ ref) ** 2)))
    cov = s / (x.size - 3) * np.linalg.inv(X.T @ X)
    assert np.allclose(se, np.sqrt(np.diag(cov)))
    # linear truth: curvature coefficient not significant
    assert abs(c[0] / se[0]) < 3
    # quadratic truth: clearly significant
    yq = 3 * x**2 + y
    q2 = fit.quadratic_fit(x, yq)
    assert abs(q2["coef"][0] / q2["se"][0]) > 10


def test_grid_ranges_from_data():
    x, y = synthetic(noise=0.05)
    a, b = fit.fit_own(x, y)
    a_lo, a_hi, b_lo, b_hi = fit.grid_ranges(a, b, x, y)
    assert a_lo < a < a_hi and b_lo < b < b_hi
    assert a_hi - a == pytest.approx(a - a_lo) and b_hi - b == pytest.approx(b - b_lo)
    # box must contain the level curve S = 4 * S_min (bounding box of that ellipse)
    Hinv = np.linalg.inv(fit.hessian(x))
    d = 3 * fit.S(a, b, x, y)
    assert a_hi - a >= math.sqrt(2 * d * Hinv[0, 0])
    assert b_hi - b >= math.sqrt(2 * d * Hinv[1, 1])
    # and not absurdly larger (at most 3x that box)
    assert a_hi - a <= 3 * math.sqrt(2 * d * Hinv[0, 0])


def test_plots_write_all_figures(tmp_path):
    x, y = synthetic(noise=0.05)
    a, b = fit.fit_own(x, y)
    fit.plots(x, y, a, b, tmp_path)
    names = sorted(p.name for p in tmp_path.glob("*.png"))
    assert names[0] == "0_dispersion.png"
    assert len(names) == 7
