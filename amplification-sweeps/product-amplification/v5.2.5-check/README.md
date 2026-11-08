# The v5.2.5 source identity

The current kernel reference is tag `v5.3.0`. This check records the source
identity at `v5.2.5`, the release the paper cites as provenance (its measured
figures are the v5.1.1 run's): every evaluated file (`include/boys` and `src`)
compared by content, one line per file, and every difference named.

The result is 73 of the 76 evaluated files byte-identical; the three that
differ are `include/boys/version.hpp` (the version constant), the
build-defaults table (which the sweeps do not read, because they name every
axis of their policy) and the coefficient header (whose only change is the
rational region-B seed, which a Chebyshev policy does not read). No file the
measurements compile changes the arithmetic they execute.

## Files

| File | What it records |
| --- | --- |
| `source-identity-v525.txt` | the tag, the counts and each differing file with its two digests and why it does not affect the measurements |
| `v5.2.5-sources-compare.txt` | the per-file comparison against v5.2.4 and v5.1.1, with each file's sha256 at all three |
| `evaluated-sources-compare-v525.txt` | the per-file comparison against v5.1.1 alone, by blob id |
| `evaluation-sources-v525.sha256` | the 76 evaluated files' sha256 at the tag |
| `version-hpp.diff` | the version header's diff, v5.1.1 to v5.2.5 |
| `make_check_525.py` | the generator that wrote them |

## Reproducing

From this directory, with the boys repository at `C:\Users\myaml\dev\_boys-library`:

    python -I make_check_525.py

It reads the tag objects and the two trees with git and rewrites the records
beside it. A tag that has moved, or a tree that differs from the tag's, fails
the comparison rather than rewriting it.
