import math

import numpy as np
import pytest

from find50.curve import crossing_log, fit_4pl, four_pl, ic50, interpolated_log_ic50

DOSES = np.array([0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0])


def test_fit_recovers_parameters_on_exact_data():
    y = four_pl(np.log10(DOSES), 5.0, 98.0, math.log10(0.4), 1.3)
    fit = fit_4pl(DOSES, y)
    assert fit.bottom == pytest.approx(5.0, abs=1e-4)
    assert fit.top == pytest.approx(98.0, abs=1e-4)
    assert fit.ec50 == pytest.approx(0.4, rel=1e-5)
    assert fit.hill == pytest.approx(1.3, rel=1e-5)
    assert fit.r2 == pytest.approx(1.0)


def test_absolute_ic50_differs_from_inflection_when_bottom_is_not_zero():
    # bottom 30 %, top 100 %: the curve crosses 50 % above the inflection point
    lec, h = math.log10(0.5), 1.0
    expected = lec + math.log10((100 - 50) / (50 - 30)) / h
    assert crossing_log(30, 100, lec, h) == pytest.approx(expected)
    res = ic50(DOSES, four_pl(np.log10(DOSES), 30, 100, lec, h))
    assert res.method == "4PL"
    assert res.ic50 == pytest.approx(10 ** expected, rel=1e-4)
    assert res.ec50_fit == pytest.approx(0.5, rel=1e-4)


def test_curve_that_never_reaches_50_has_no_crossing():
    assert math.isnan(crossing_log(60, 100, 0.0, 1.0))


def test_interpolated_crossing():
    conc = [0.1, 1.0, 10.0]
    li, q = interpolated_log_ic50(conc, [90, 70, 30])
    assert q == "" and li == pytest.approx(0.0 + (70 - 50) / (70 - 30))
    assert interpolated_log_ic50(conc, [40, 30, 20]) == (pytest.approx(-1.0), "<")
    assert interpolated_log_ic50(conc, [95, 80, 60]) == (pytest.approx(1.0), ">")


def test_out_of_range_results_are_qualified():
    low = ic50(DOSES, four_pl(np.log10(DOSES), 0, 100, math.log10(0.001), 1))
    high = ic50(DOSES, four_pl(np.log10(DOSES), 0, 100, math.log10(1000), 1))
    assert (low.qualifier, low.ic50, low.label) == ("<", 0.01, "< 0.01")
    assert (high.qualifier, high.ic50, high.label) == (">", 10.0, "> 10")
    assert high.ic50_fit > 10                      # the extrapolated fit is still reported


def test_too_few_points_fall_back_to_interpolation():
    res = ic50([1, 5, 25, 125], [80, 55, 25, 22])
    assert res.method == "interpolation"
    assert any("too few" in f for f in res.flags)


def test_poor_fit_falls_back_to_interpolation():
    y = np.array([100, 40, 95, 30, 90, 20, 10.0])
    res = ic50(DOSES, y)
    assert res.method == "interpolation"
    assert any("poor 4PL fit" in f for f in res.flags)
    assert any("non-monotonic" in f for f in res.flags)


def test_noisy_data_ci_contains_truth():
    rng = np.random.default_rng(3)
    y = four_pl(np.log10(DOSES), 0, 100, math.log10(0.2), 1.0) + rng.normal(0, 2, DOSES.size)
    res = ic50(DOSES, y)
    assert res.method == "4PL"
    assert res.ci95_low < 0.2 < res.ci95_high
    assert res.ic50 == pytest.approx(0.2, rel=0.2)


def test_bad_input():
    assert ic50([1.0], [50.0]).method == "none"
    assert ic50([0, 0, 1], [100, 90, 20]).n_points == 1     # zero doses are not doses
