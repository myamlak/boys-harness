# Materiality saturation

`b7_ratio_distribution.json` reduces the deposited v5.0.0 cube to the 1,808
four-cell contrasts of the resolving lanes and reports the distribution of
each contrast's largest signed interaction-to-response ratio. The ratios are
not near a threshold: they are saturated above it.

On the fused route, 828 of the 1,648 informative contrasts have a ratio above
1, and the median is `1.0010077181740726`; the upper quantiles are
`1.966900515945531` (75 per cent), `1.9996621529640912` (90 per cent) and
`3.9264180858800426` (maximum). On the separate route, 836 are above 1 and
the median is `1.037219146449483`. Both routes have the same two neighbouring
bands: 11 contrasts in `(0.05, 0.10]` and 2 in `(0.10, 0.20]`.

That is why the material counts barely move with the margin. Raising it from
5 to 10 per cent removes 11 contrasts (fused `912 -> 901`, separate
`914 -> 903`), and raising it from 10 to 20 per cent removes 2 more (fused
`901 -> 899`, separate `903 -> 901`). The 5 per cent margin already leaves
more than half the informative contrasts above 1, so the audit's materiality
verdict does not rest on the choice of margin.

The saturation survives the v5.3.3 re-measurement: the changed fp64 cells
lower the fused counts to `897 / 885 / 883` and the separate counts to
`895 / 891 / 885`, with the same verdict at every margin. The absolute-value
(cell-max) reproduction in the same reduction is a different statistic and is
much lower (`466 / 388 / 332` fused, `474 / 418 / 374` separate); this note
concerns the signed pointwise statistic the audit uses.
