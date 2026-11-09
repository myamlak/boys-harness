# C strong optional

Two independent measurements, deposited with their scripts, records and generated
report:

- `best-vs-best/` - a best-uniform (minimax) polynomial arm against the corrected
  BRASIL rational record, under matched interval, target, 70-digit reference,
  binary64 monomial storage and fused/separate Horner evaluation on the 2,001-
  and 4,001-node grids. Stored counts are reported separately from the paper's
  existing pipeline comparison.
- `root-isolation/` - exact Sturm root isolation and certified interval-wide
  denominator bounds for all 98 accepted rational fits in the corrected x0,
  second-edge a2 and 60-digit control records.

`REPORT.md` is generated from the two records by `summarize.py`; every number in
it is read from the records. Selftests:

```sh
python root-isolation/root_isolation.py --selftest
python best-vs-best/best_vs_best.py --selftest
python summarize.py
```
