"""Find50 reproduces the IC50 definition of the original notebook.

`legacy_ic50` is a port of `GetData` / `FittingMath` / `Grapher` from
legacy/IC50 FINAL PROGRAM-5.ipynb: normalise to the maximum (= 100), drop rising
points with the notebook's fixed-index rules, fit a Hill-1 logistic with
`scipy.optimize.leastsq` (for R^2), then step s = log10(c) up from -2 (240 steps of
0.01, then 84 of 0.04) until the linearly interpolated data fall to <= 50. The port is
checked against output recorded in the original notebooks before comparing."""
import numpy as np
import pytest
from scipy.optimize import leastsq

from find50.curve import fit_4pl, interpolated_log_ic50

LEGACY_X = [-2, -1, -0.69897000433, -0.39794000867, -0.096910013, 0.20411998265, 0.50514997832]


def normalize(arr):
    return [(i / max(arr)) * 100 for i in arr]


def legacy_clean(y0):
    """GetData's point handling (returns x, y used for the fit and the crossing)."""
    x = list(LEGACY_X)
    y = normalize(y0)
    if y[6] > y[5]:
        y.pop(6), x.pop(6)
        y = normalize(y)
    if y[5] > y[4]:
        y.pop(5), x.pop(5)
        y = normalize(y)
    if y[4] > y[3]:
        y.pop(3), x.pop(3)
    if y[3] > y[2]:
        y.pop(3), x.pop(3)
        y = normalize(y)
    if y[2] > y[1]:
        y.pop(2), x.pop(2)
    if y[1] > y[0]:
        y[1] = y[0]
    return x, y


def legacy_fit_r2(x, y):
    eq = lambda v, x: v[0] + (v[1] - v[0]) / (1 + 10 ** (v[2] - np.asarray(x)))  # noqa: E731
    err = lambda v, x, y: eq(v, x) - np.asarray(y)                                # noqa: E731
    out = leastsq(err, [np.min(y), np.max(y), np.mean(x)], args=(x, y), full_output=1)
    fvec = out[2]["fvec"]
    return 1 - np.sum(fvec ** 2) / np.sum((np.asarray(y) - np.mean(y)) ** 2)


def legacy_crossing(x, y, steps=((0.01, 240), (0.04, 84))):
    s = -2
    g = np.interp(s, x, y)
    for step, n in steps:
        for _ in range(n):
            if g > 50:
                s += step
                g = np.interp(s, x, y)
    return s


def legacy_ic50(y0):
    x, y = legacy_clean(y0)
    if min(y) > 50:
        return ">>3.2", None
    return 10 ** legacy_crossing(x, y), legacy_fit_r2(x, y)


# responses and printed results recorded in the original notebook (first two compounds)
RECORDED = [
    ([0.23692306, 0.27338662, 0.10430791, 0.07424567, 0.0794591, 0.06991952, 0.06709346],
     0.16982436524617486, 0.8521481176821908),
    ([0.26199642, 0.26588794, 0.05980483, 0.03407334, 0.03364215, 0.02345945, 0.02129781],
     0.15848931924611173, 0.8175698991804989),
    ([0.2372854, 0.25468987, 0.08843646, 0.07850028, 0.06818543, 0.06675838, 0.06967394],
     0.16982436524617486, 0.8064780423506968),
    ([0.25921782, 0.27413032, 0.03921582, 0.02746651, 0.02109272, 0.01916445, 0.01848225],
     0.1479108388168211, 0.7928547557474162),
]


@pytest.mark.parametrize("y0, ic50_printed, r2_printed", RECORDED)
def test_port_reproduces_recorded_notebook_output(y0, ic50_printed, r2_printed):
    v, r2 = legacy_ic50(y0)
    assert v == pytest.approx(ic50_printed, rel=1e-9)
    assert r2 == pytest.approx(r2_printed, rel=1e-6)


@pytest.mark.parametrize("y0, ic50_printed, r2_printed", RECORDED)
def test_same_crossing_and_fit_on_the_same_points(y0, ic50_printed, r2_printed):
    x, y = legacy_clean(y0)
    li, q = interpolated_log_ic50(10 ** np.asarray(x), y)
    # the original stops at the first 0.01 step at or past the exact crossing
    assert q == "" and 0 <= np.log10(ic50_printed) - li < 0.01 + 1e-9
    fit = fit_4pl(10 ** np.asarray(x), y, fixed_hill=1.0)
    assert fit.r2 == pytest.approx(r2_printed, abs=2e-3)


@pytest.mark.parametrize("seed", range(20))
def test_random_curves_agree_within_one_legacy_step(seed):
    rng = np.random.default_rng(seed)
    lec = rng.uniform(-0.9, 0.4)
    x = np.array(LEGACY_X[1:])
    y = 100 / (1 + 10 ** (rng.uniform(0.7, 2) * (x - lec))) + rng.normal(0, 3, x.size)
    y = np.maximum.accumulate(y[::-1])[::-1]               # monotonic, as GetData enforces
    xs, ys = [-2] + list(x), [max(100, y[0])] + list(y)
    s = legacy_crossing(xs, ys)
    li, q = interpolated_log_ic50(10 ** np.asarray(xs), ys)
    if q == "":
        step = 0.01 if s <= 0.4 + 1e-9 else 0.04
        assert 0 <= s - li < step + 1e-9
