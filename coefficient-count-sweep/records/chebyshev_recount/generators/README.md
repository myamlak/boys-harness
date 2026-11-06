# The two coefficient generators

These are the two revisions of the `boys` repository's
`tools/gen_boys_coefficients.py` that the recount record names. Each file is
the tag's blob, byte for byte.

| File | `boys` tag | Commit | SHA-256 |
| --- | --- | --- | --- |
| `gen_boys_coefficients_v5.0.0.py` | `v5.0.0` | `482fd3388c276091834a6c12a642d9de2481c9f2` | `ccc8b011f4170ce0b530bc4781b09670f42839fa18caaa38ea7d3a9ea9ba9bb4` |
| `gen_boys_coefficients_v5.1.1.py` | `v5.1.1` | `834bdbb18a30ac8538b519436429f78d41f187ce` | `11185fba03abe74d01b30f093818f48972eb8d88413d446b505198b869115180` |

The earlier counts came from the v5.0.0 generator
(`boys.previous_generator_sha256` in `records/chebyshev_recount.json`). The
corrected recount re-solved the Chebyshev side with the v5.1.1 generator
(`boys.generator_sha256`), the revision that corrects the split-Clenshaw odd
seed.

The v5.2.5 generator is a different source: SHA-256
`1d8b88ac30c757138a4c87219d68b56756c175c36230818b2360dc0bf7bdd373`. It is
the generator named by the records re-measured at v5.2.5; it is neither
generator the recount record names and is not deposited here.

The recount record gives `boys.version` 5.1.1 and `describe`
`v5.1.0-1-g3b187d5`. That describe hash is not a commit in the deposited
`boys` history. The generator digest is the identity the record uses, and the
tagged v5.1.1 file has that digest.
