# Per-call error over candidate per-axis sums

`b1_calls.tsv` reduces the deposited v5.0.0 cube to one row per joint call:
576 calls per configuration, 1,152 rows. `max_abs_err` is the call's largest
`abs_err` over its 858 points. The ratio columns divide that error by four
candidate denominators:

- `fig_ratio`: the accessor's figure for the call's cell and point;
- `statement_ratio`: the row's whole published statement, which adds the
  value-dependent half digit on the fp16 and bf16 lanes;
- `six_statement_ratio`: six copies of that statement;
- `six_base_ratio` and `base_ratio`: the same over the lane base alone,
  subtracting the 1e-7 form term on the float-like lanes.

## Maxima

| configuration | lane | max err / figure | max err / statement | max err / (6 x statement) |
| --- | --- | ---: | ---: | ---: |
| fused | fp64 | 0.629268025193091 | 0.629268025193091 | 0.1048780041988485 |
| fused | fp32 | 0.8191245254432976 | 0.8191245254432976 | 0.1365207542405496 |
| fused | fp16 | 1614.0104161801694 | 0.9934027438291815 | 0.16556712397153023 |
| fused | bf16 | 8306.411848962533 | 0.9985908700338076 | 0.16643181167230128 |
| separate | fp64 | 0.6270981920683498 | 0.6270981920683498 | 0.10451636534472497 |
| separate | fp32 | 0.8191245254432976 | 0.8191245254432976 | 0.1365207542405496 |
| separate | fp16 | 1614.0104161801694 | 0.9934027438291815 | 0.16556712397153023 |
| separate | bf16 | 8306.411848962533 | 0.9985908700338076 | 0.16643181167230128 |

The largest ratio anywhere is `0.9985908700338076` over the whole statement at
bf16, and `0.8191245254432976` over the figure on the resolving lanes. On the
half lanes the accessor figure alone is far smaller than the delivered error,
because the row's statement adds the half digit the accessor figure does not
carry; their figure ratios are therefore not bounds and are listed only to
show why the statement is the right denominator.

## What the sum is not

The exact sum over six separate per-axis figures is not defined by the
published rows. The accessor's figure already folds the lane base and the
conditional terms into one number, while the six published axis rows are
heterogeneous: some state the lane base, some the accessor figure, and the
half lanes add a value-dependent half digit that no separate axis row
publishes. Six copies of a row are therefore not six per-axis figures, and a
quotient against them is a deliberately generous surrogate, not a measured
sum.

What the record supports is the table above: on the resolving lanes the
measured call error stays below the accessor figure (`0.819` at worst), and
on every lane it stays below the row's whole statement (`0.999` at worst).
Even six whole statements leave at least a factor of six in hand
(`0.1665` at worst).
