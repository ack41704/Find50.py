# Original notebook

`IC50 FINAL PROGRAM-5.ipynb` is the original Find50 notebook. Its outputs have been
cleared; the code is unchanged. Keep it to reproduce earlier results exactly. For new
analyses, use the `find50` package. [`tests/test_legacy_equivalence.py`](../tests/test_legacy_equivalence.py)
contains a port of this notebook's IC50 calculation, which reproduces its recorded
output.

## What it does

- **Cell 1** defines `GetData`, `FittingMath`, `residuals` and `Grapher`. The response
  series is normalised to its maximum, and rising points are removed with fixed-index
  rules. A logistic with a Hill slope of 1 is fitted (`scipy.optimize.leastsq`); the
  fit is used for R² only. The IC50 is then found by stepping log10(c) up from −2 in
  0.01 steps (240 steps, then 84 steps of 0.04) until the linearly interpolated data
  fall to 50 %. It prints `10**s` and R², or `>>3.2` when no point is below 50 %.
- **Cell 2** reads `2023.71.xlsx` with `pandas.read_excel`. It takes the association
  readout from data row 154 (75 s) and the dissociation readout from row 264 (130 s),
  and loops over 8 compounds in 16-column blocks. Within each block, column 1 is used as the
  no-inhibitor point (placed at 0.01 µM), column 3 is the positive-control baseline,
  and columns 5–15 are 0.1–3.2 µM. "Rel" variants subtract the smaller of the control and
  highest-dose responses.

## Known issues (fixed in the package)

- The file name, row indices, block stride, number of compounds and concentrations are
  hard-coded. Files with fewer than 8 compounds end with an `IndexError`.
- In `GetData`, a rise at index 4 removes index 3. The no-inhibitor well is treated as
  a dose at 0.01 µM, so a curve that is already below 50 % at the first dose prints
  `0.01`.
- In the dissociation "Rel" branch, the highest-dose baseline is read from the
  association row (`bd2 = a[15 + e]`). When the control is lower, the highest dose is
  set to 0. If the two baselines are equal, the absolute series is used unchanged.
- The IC50 comes from the step search on the data, not from the fitted curve. Its
  resolution is 0.01 log units (2.3 %), and 0.04 log units above about 2.5 µM.
- `shapely` is imported but not used.
