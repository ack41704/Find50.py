"""
Concentration-response curves and IC50 values.

Responses are expressed as **% activity**: 100 % = no inhibitor (top reference), 0 % =
no signal (absolute) or the positive-control level (relative). Two IC50 estimates are
computed for every curve:

* ``ic50_fit``    — the concentration at which the fitted four-parameter logistic (4PL)
  curve crosses 50 % activity (the *absolute* IC50):

      y(x) = bottom + (top - bottom) / (1 + 10**(hill * (x - log_ec50))),   x = log10(c)

      log_ic50 = log_ec50 + log10((top - 50) / (50 - bottom)) / hill

  ``ec50_fit`` = 10**log_ec50 is the inflection point of the curve.
* ``ic50_interp`` — the model-free crossing: the lowest concentration at which the
  measured points, linearly interpolated in log10 concentration, fall to 50 %. This is
  the definition used by the original notebook.

The reported ``ic50`` is the 4PL value when the fit is acceptable (R^2 >= r2_min and the
crossing lies inside the tested range), otherwise the interpolated value. If every dose
is below 50 % the result is ``< lowest concentration``; if none is, ``> highest
concentration``.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
import warnings

import numpy as np


def four_pl(x, bottom, top, log_ec50, hill):
    """4PL in log10 concentration (decreasing for hill > 0)."""
    return bottom + (top - bottom) / (1.0 + 10.0 ** (hill * (np.asarray(x, float) - log_ec50)))


def crossing_log(bottom, top, log_ec50, hill, level=50.0):
    """log10 concentration where the 4PL equals `level` (nan if it never does)."""
    lo, hi = min(bottom, top), max(bottom, top)
    if not lo < level < hi or hill == 0:
        return math.nan
    return log_ec50 + math.log10((top - level) / (level - bottom)) / hill


@dataclass
class FourPLFit:
    bottom: float
    top: float
    log_ec50: float
    hill: float
    r2: float
    log_ec50_se: float = math.nan

    @property
    def ec50(self):
        return 10.0 ** self.log_ec50

    def log_ic50(self, level=50.0):
        return crossing_log(self.bottom, self.top, self.log_ec50, self.hill, level)

    def __call__(self, x):
        return four_pl(x, self.bottom, self.top, self.log_ec50, self.hill)


def fit_4pl(conc, activity, hill_bounds=(0.1, 10.0), fixed_hill=None) -> FourPLFit | None:
    """Least-squares 4PL fit of % activity vs concentration (concentrations > 0).

    Bounds keep the fit physical: top in [0, 200] %, bottom in [-50, 100] %, Hill slope
    in `hill_bounds` (positive = activity decreases with dose), log EC50 within 2 log
    units of the tested range. Returns None with fewer than 4 points."""
    from scipy.optimize import curve_fit

    c = np.asarray(conc, float)
    y = np.asarray(activity, float)
    ok = np.isfinite(c) & np.isfinite(y) & (c > 0)
    x, y = np.log10(c[ok]), y[ok]
    n_par = 3 if fixed_hill is not None else 4
    if len(x) < n_par or np.ptp(x) == 0:
        return None
    xlo, xhi = x.min() - 2, x.max() + 2
    p0 = [float(np.clip(y.min(), -49, 99)), float(np.clip(y.max(), 1, 199)),
          float(np.median(x)), 1.0 if fixed_hill is None else fixed_hill]
    lower = [-50.0, 0.0, xlo, hill_bounds[0]]
    upper = [100.0, 200.0, xhi, hill_bounds[1]]
    # start the inflection where the data cross half-way between their extremes
    mid = (y.max() + y.min()) / 2
    order = np.argsort(x)
    below = np.nonzero(y[order] <= mid)[0]
    if len(below):
        p0[2] = float(np.clip(x[order][below[0]], xlo, xhi))
    if fixed_hill is None:
        f, start, lo, hi = four_pl, p0, lower, upper
    else:
        f = lambda xx, b, t, e: four_pl(xx, b, t, e, fixed_hill)  # noqa: E731
        start, lo, hi = p0[:3], lower[:3], upper[:3]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            popt, pcov = curve_fit(f, x, y, p0=start, bounds=(lo, hi), maxfev=20000)
    except (RuntimeError, ValueError):
        return None
    pred = f(x, *popt)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan
    se = math.nan
    if len(x) > n_par and np.all(np.isfinite(pcov)):
        se = float(np.sqrt(max(pcov[2, 2], 0.0)))
    hill = float(popt[3]) if fixed_hill is None else float(fixed_hill)
    return FourPLFit(float(popt[0]), float(popt[1]), float(popt[2]), hill, r2, se)


def interpolated_log_ic50(conc, activity, level=50.0):
    """(log10 crossing, qualifier) from linear interpolation of the points in log10
    concentration, scanning from the lowest dose: ("<" when the lowest dose is already
    at or below `level`, ">" when no dose reaches it)."""
    c = np.asarray(conc, float)
    y = np.asarray(activity, float)
    ok = np.isfinite(c) & np.isfinite(y) & (c > 0)
    x, y = np.log10(c[ok]), y[ok]
    if not len(x):
        return math.nan, ""
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    if y[0] <= level:
        return float(x[0]), "<"
    for i in range(1, len(x)):
        if y[i] <= level:
            if y[i] == y[i - 1]:
                return float(x[i]), ""
            return float(x[i - 1] + (level - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1])), ""
    return float(x[-1]), ">"


@dataclass
class IC50Result:
    ic50: float                 # reported IC50 (same unit as the concentrations); see qualifier
    qualifier: str              # "" | "<" | ">"
    method: str                 # "4PL" | "interpolation" | "range"
    ic50_fit: float = math.nan  # absolute IC50 of the 4PL fit
    ec50_fit: float = math.nan  # inflection point of the 4PL fit
    ic50_interp: float = math.nan
    ci95_low: float = math.nan  # 95 % CI of the 4PL inflection (log-normal, from the fit SE)
    ci95_high: float = math.nan
    hill: float = math.nan
    top: float = math.nan
    bottom: float = math.nan
    r2: float = math.nan
    n_points: int = 0
    flags: list = field(default_factory=list)
    fit: FourPLFit | None = field(default=None, repr=False)

    @property
    def label(self) -> str:
        """Text for tables: "0.861", "< 0.1", "> 3.2" (3 significant figures)."""
        if not np.isfinite(self.ic50):
            return "n/a"
        v = f"{self.ic50:.3g}"
        return f"{self.qualifier} {v}" if self.qualifier else v

    def as_dict(self) -> dict:
        d = asdict(self)
        d.pop("fit")
        d["flags"] = "; ".join(self.flags)
        d["label"] = self.label
        return d


def ic50(conc, activity, r2_min: float = 0.8, level: float = 50.0,
         fixed_hill: float | None = None, min_points_fit: int = 5) -> IC50Result:
    """IC50 of one concentration series (% activity, concentrations > 0).

    The 4PL value is reported only with at least `min_points_fit` doses (a free 4PL
    through 4 points always fits exactly), R^2 >= `r2_min` and a 50 % crossing inside
    the tested range; otherwise the interpolated crossing is reported."""
    c = np.asarray(conc, float)
    y = np.asarray(activity, float)
    ok = np.isfinite(c) & np.isfinite(y) & (c > 0)
    c, y = c[ok], y[ok]
    flags = []
    if len(c) < 2:
        return IC50Result(math.nan, "", "none", n_points=len(c), flags=["fewer than 2 points"])
    cmin, cmax = float(c.min()), float(c.max())
    order = np.argsort(c)
    if np.any(np.diff(y[order]) > 10.0):
        flags.append("non-monotonic (activity rises > 10 points between doses)")

    li, q = interpolated_log_ic50(c, y, level)
    interp = 10.0 ** li if np.isfinite(li) else math.nan
    fit = fit_4pl(c, y, fixed_hill=fixed_hill)
    res = dict(ic50_interp=interp if not q else math.nan, n_points=len(c), fit=fit)
    if fit is not None:
        lic = fit.log_ic50(level)
        res.update(ic50_fit=10.0 ** lic if np.isfinite(lic) else math.nan, ec50_fit=fit.ec50,
                   hill=fit.hill, top=fit.top, bottom=fit.bottom, r2=fit.r2)
        if np.isfinite(fit.log_ec50_se):
            res.update(ci95_low=10.0 ** (fit.log_ec50 - 1.96 * fit.log_ec50_se),
                       ci95_high=10.0 ** (fit.log_ec50 + 1.96 * fit.log_ec50_se))
        enough = len(c) >= (min_points_fit if fixed_hill is None else min_points_fit - 1)
        good = np.isfinite(fit.r2) and fit.r2 >= r2_min
        in_range = np.isfinite(lic) and cmin <= 10.0 ** lic <= cmax
        if not enough:
            flags.append(f"{len(c)} doses: too few to report a 4PL IC50")
        elif not good:
            flags.append(f"poor 4PL fit (R2 = {fit.r2:.2f} < {r2_min})")
        if enough and good and in_range and not q:
            return IC50Result(10.0 ** lic, "", "4PL", flags=flags, **res)
    else:
        flags.append("4PL fit not possible (fewer than 4 points or no convergence)")
    if q == "<":
        return IC50Result(cmin, "<", "range", flags=flags + [f"all doses <= {level:g} %"], **res)
    if q == ">":
        return IC50Result(cmax, ">", "range", flags=flags + [f"no dose reaches {level:g} %"], **res)
    return IC50Result(interp, "", "interpolation", flags=flags, **res)
