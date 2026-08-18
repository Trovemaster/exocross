# ExoCross benchmark comparator

`benchmark_exocross.py` reads an ExoCross input, determines the expected
cross-section outputs, and compares newly generated files with stored reference
files. It uses only the Python standard library.

The first supported task is `absorption` with `voigt-super`, including
`temperature-list` and `pressure-list`. For every pressure, the comparator
expects the ExoCross filename

```text
<output>__<pressure-as-1.00000000E+000>.xsec
```

and validates that every row has one wavenumber plus one cross-section column
per temperature, in input-list order.

## Run

Generate fresh outputs in a directory separate from the stored references,
then run:

```powershell
python benchmark_exocross.py `
  C:\sergei\programs\exocross-inputs\benchmarks\CO_01\CO_T-list_P-list_V_01.inp `
  --candidate-dir C:\path\to\fresh\CO_01
```

The reference directory defaults to the input file's directory. Use
`--reference-dir` to override it. The process exits with code 0 for a match, 1
for numerical/shape mismatches, and 2 for invalid input or unreadable files.
`--json` provides output suitable for CI.

By default, cross sections must agree to six significant figures. This is
implemented as a symmetric relative tolerance of `5e-6`, with no ordinary
absolute tolerance, so a discrepancy around `1e-30` is not hidden merely
because both values are small. A `1e-300` relative scale floor only makes zero
and subnormal handling well-defined. The wavenumber grid is checked separately
to nine significant figures.

Useful controls are:

```text
--significant-figures N
--grid-significant-figures N
--absolute-tolerance VALUE
--relative-floor VALUE
--max-examples N
```

An absolute tolerance should be enabled only when the physical benchmark has a
justified noise floor.

## Tests

```powershell
python -m unittest -v
```
