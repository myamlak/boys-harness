# Neutral copies

The records in this deposit name an internal tree of a sibling manuscript in
one metadata field. These copies present the same records without that name.
No measured value, table row or digest is changed; the only fields changed
are the ones listed below.

| Copy | Original | Change |
| --- | --- | --- |
| `fullcube/fullcube_decomposition.py` | `composition-audit/fullcube/fullcube_decomposition.py` | the docstring names the full-cube decomposition lane instead of the sibling manuscript |
| `signed-interaction/signed_interaction.py` | `composition-audit/signed_interaction.py` | the emitted `paths_relative_to` is `the original run root`; the import search path also names this repository's `composition-audit` directory |
| `signed-interaction/signed-interaction.json` | `composition-audit/records/signed-interaction.json` | `paths_relative_to` is `the original run root` |
| `synthetic-surface/synthetic_surface.py` | `composition-audit/synthetic_surface.py` | the emitted `paths_relative_to` is `the original run root` |
| `synthetic-surface/synthetic-surface.json` | `composition-audit/records/synthetic-surface.json` | `paths_relative_to` is `the original run root` |

In the copies, `paths_relative_to` is a description, not a path. The other
path fields are unchanged: they name the tree the measurement ran in, whose
cell files are not deposited. The deposited equivalents are
`composition-audit/records/cells/` and `composition-audit/records/reference.csv`.

The originals were not edited in place because deposited records bind their
bytes:

- `composition-audit/fullcube/fullcube_summary.json` records the sha256 of
  `fullcube_decomposition.py`.
- `composition-audit/records/synthetic-surface.json` records the sha256 of the
  deposited `synthetic_surface.py`.
- `composition-audit/records/figure-ratio.json` records the sha256 of
  `signed-interaction.json` and documents that the deposited
  `signed_interaction.py` regenerates that record.
- `paths_relative_to` is the base for every path field in the two records, so
  a value that does not describe the run tree would make the record misstate
  where its inputs live.

Use these copies when the internal name must not appear. The originals remain
the authority for the measured values.
