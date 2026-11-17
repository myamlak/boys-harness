# Neutral copies

The campaign studies cite the copies in this directory. Each copy preserves
the relative path of its original below `composition-audit/`.

The data and digests are identical. The copy is the citation target. The
original record carries a historical internal name.

| Copy | Original |
| --- | --- |
| `neutral/compare_stability.py` | `compare_stability.py` |
| `neutral/fullcube/fullcube_decomposition.py` | `fullcube/fullcube_decomposition.py` |
| `neutral/gen_reference.py` | `gen_reference.py` |
| `neutral/records/signed-interaction.json` | `records/signed-interaction.json` |
| `neutral/records/synthetic-surface.json` | `records/synthetic-surface.json` |
| `neutral/selftest_stability.py` | `selftest_stability.py` |
| `neutral/selftest_summarise.py` | `selftest_summarise.py` |
| `neutral/signed_interaction.py` | `signed_interaction.py` |
| `neutral/summarise_interaction.py` | `summarise_interaction.py` |
| `neutral/synthetic_surface.py` | `synthetic_surface.py` |

The copies carried from v1.0.41 under `neutral/signed-interaction/` and
`neutral/synthetic-surface/` are consolidated at the mirror paths here;
the v1.0.41 tag retains the earlier paths.

These are citation artifacts. Do not run `neutral/synthetic_surface.py` in
place: it resolves `records/` relative to `composition-audit/` and would
overwrite the deposited record.
