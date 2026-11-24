# Clean-edge solver independence

This directory records an independent re-solve of the two clean-edge cells
whose 60-digit counts carry the deep-target result: the `2^5` and `ship`
widths at edge `a = 62` and relative target `1e-15`. The question is whether
the 60-digit meeting counts, 15 at `2^5` and 13 at `ship`, are specific to
BRASIL/`baryrat`. A from-scratch rational Remez/Newton solver and an iterated
weighted-least-squares (Sanathanan-Koerner) construction, both at 70 decimal
digits, reach the same counts; de la Vallee Poussin bounds from independent
fits exclude the count one below for every admissible pair.

The counts under test are the control's, not the binary64 record's. The
record stores 22 at `2^5` and 16 at `ship`; the 60-digit control deposited in
`../clean-edge-60dig-2026-11-24/` lowers them to 15 and 13. This deposit tests
those two numbers.

## Protocol

The record is `../clean-edge-2026-11-22/edge62_sweep.json` (SHA-256
`963facd9e18c96b547b6543c8713d5dd8294104d2776fb8d71dd30c0d40bff8b`). The
two accepted control cells are
`../clean-edge-60dig-2026-11-24/cells/p5_rel_1e-15.json` (SHA-256
`b7dd16a6172237b250092dc6c52c0708b4514b5af790eb04bf5f821cb60ba052`) and
`../clean-edge-60dig-2026-11-24/cells/ship_rel_1e-15.json` (SHA-256
`70cf661f15f9f4cbc9d5665b2df40cebb65c549bf8ae7905809ade6cc845455b`). The
same digests are carried in `out/summary.json`.

- Edge `a = 62.0`; record widths `W = 166.72481653881133` (`2^5`) and
  `W = 99.51466602197347` (`ship`); `t` uniform on `[-1, 1]`,
  `x = a + W (t + 1)/2`.
- Reference `F_0(x) = sqrt(pi)/(2 sqrt(x)) erf(sqrt(x))`, evaluated at 70
  decimal digits.
- 2,001-node certifying grid (2,000 intervals), built in binary64 as the
  control builds it; a 4,001-node check grid; a 20,001-node grid for the
  final fits.
- Relative target `1e-15`: `|P/Q - F_0|/F_0 <= 1e-15` at every node.
- Rational type `[m/k]` with `m >= 1` and denominator degree `1 <= k <= 12`;
  stored count `m + k + 1`; every pair at a stored count is swept.
- Denominator guard: `min |Q| / max |Q| >= 1e-3` on the grid, the control's
  dimensionless guard.
- Coefficients held at 60 significant digits and evaluated at 60 decimal
  digits; the solver, the coefficient conversion and `F_0` use 70 decimal
  digits.
- A stored count is the lowest count at which an admissible pair meets.
  `e60` is the maximum 60-digit relative error on the 2,001-node grid;
  `guard60` is the minimum-to-maximum denominator ratio on that grid.

## Independent solver

`remez_rational.py` is a from-scratch classical rational Remez exchange at 70
decimal digits with `mpmath`. It uses differential correction (30 iterations)
followed by an exact Newton polish (40 iterations) on the reference
equations; reference points are continuous, refined by parabolic
interpolation, and each exchange keeps the `N = m + k + 2` alternating extrema
with the largest minimum magnitude (global de la Vallee Poussin and Maehly
exchanges). It declares convergence only on the alternation-theorem
certificate: `N` alternating reference extrema of common magnitude with the
grid maximum within `1e-10` relative. The weight is `1/F_0`, so the minimised
quantity is exactly the control's relative error. The solver is independent
of BRASIL/`baryrat` (barycentric Lawson iteration) and of the record's
binary64 kernel.

Known-answer validation (`--test`): the minimax polynomial of degree `n-1`
for `t^n` on `[-1, 1]` is `t^n - T_n/2^(n-1)`, with error `2^(1-n)`. The
solver recovers `t^4`, `t^5`, `t^6` and `t^8` with those errors to `1e-12`
relative, and the `t^4` coefficients equal `t^2 - 1/8` to `1e-40`. The test
passes in 1.0 s.

No suitable rational-minimax package was found on PyPI; `chebfun` 0.11.0 is
the polynomial `chebpy` package, and the other candidate names are not
rational minimax solvers. The independent arm is therefore the from-scratch
implementation.

The constructive arm (`sk_pairs.py`) is the iterated weighted-least-squares
(Sanathanan-Koerner) fit at 70 decimal digits, 15 iterations per pair. It
supplies admissible fits at every pair and initialises the Newton/Remez
polish. A cold Remez solve of `2^5` `[6/8]` did not certify: differential
correction from Chebyshev references fell into a spurious near-pole local
solution at `1.753192978686451e-15` (grid/reference ratio 131.0), and
denominator-degree continuation `[6/0]` to `[6/8]` reached the same basin. The
weighted-least-squares fit is what makes the targeted Newton start well enough
to meet.

## Results

Both counts are reproduced. The meeting fits:

| Width | Control count | Independent count | Pair | e60 at 2,001 nodes | e60 at 4,001 nodes | e60 at 20,001 nodes | guard60 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `2^5` | 15 | 15 | `[6/8]` | `2.5634851571260719e-16` | `2.5767742029803888e-16` | `2.5903440950972191e-16` | `7.8417105153113476e-3` |
| `ship` | 13 | 13 | `[6/6]` | `5.3885322603847251e-16` | `5.4026100869330508e-16` | `5.5277638127626663e-16` | `4.452579531066782e-2` |

All three grids meet at both widths; `out/final_checks.json` holds the values.
The `2^5` fit is the targeted Newton polish of the weighted-least-squares
`[6/8]` fit, whose own 60-digit score is `8.0230334461484095e-16` (guard
`7.9104785992465589e-3`) and which already meets on the 2,001- and
4,001-node grids. At `2^5` the independent arm also meets at `[7/7]`, with
`e60 = 7.0226416642394374e-16` and guard `7.5447081937891473e-3`. At `ship`
the weighted-least-squares `[6/6]` fit is `1.8067095438838784e-15` and does
not meet; the targeted Newton polish reaches `5.3885e-16`.

The control's own accepted fits were re-scored independently from their
stored coefficients (`verify_scoring.py`). The stored `e60` values reproduce
exactly at 2,001 nodes: `6.4804911767675102e-16` at `2^5` `[6/8]` (guard
`7.6084740911630768e-3`) and `8.8616641776340029e-16` at `ship` `[6/6]`
(guard `4.4287004445073072e-2`). Both stored 60-digit coefficient sets are
the 60-digit rounding of the stored 70-digit sets, and both fits meet on the
4,001- and 20,001-node grids. This is the building-block check that the
scoring route is faithful; it does not use the independent fits.

At one count below, the sweep finds no meeting pair:

| Width | Stored count | Pairs swept | Best independent e60 | Meeting pairs |
| --- | --- | --- | --- | --- |
| `2^5` | 14 | 12 | `8.4753572430358747e-15` at `[6/7]` | none |
| `ship` | 12 | 10 | `2.9284113044804476e-14` at `[5/6]` | none |

## Lower bounds

The absence of a meeting fit at the count below is not by itself an
exclusion. `lower_bounds.py` turns each independent fit into a de la Vallee
Poussin lower bound. For a fit `r` of type `[m/k]`, let `N = m + k + 2`. If
the weighted error of `r` has `N` alternating grid-node extrema of magnitude
at least `L`, then every guard-admissible rational of the same type has grid
maximum at least `L`: the difference of two such fits is `R/(Q Q')` with
`deg R <= m + k`, and alternation forces `m + k + 1` zeros of `R`. The bound
is stated on grid nodes, the conservative form for the control's grid
maximum.

| Width | Count below | Pairs | Minimum grid bound `L` | Binding pair | `L / 1e-15` | Pairs above target |
| --- | --- | --- | --- | --- | --- | --- |
| `2^5` | 14 | 12 | `1.3249223056307639e-15` | `[6/7]` | 1.3249 | 12 / 12 |
| `ship` | 12 | 10 | `5.359156349772307e-15` | `[5/6]` | 5.3592 | 10 / 10 |

The three tightest pairs at each count:

| Width | Pair | Grid bound `L` | `L / 1e-15` |
| --- | --- | --- | --- |
| `2^5` | `[6/7]` | `1.3249223056307639e-15` | 1.3249 |
| `2^5` | `[7/6]` | `1.5356701247475752e-15` | 1.5357 |
| `2^5` | `[5/8]` | `2.0434674466993266e-15` | 2.0435 |
| `ship` | `[5/6]` | `5.359156349772307e-15` | 5.3592 |
| `ship` | `[6/5]` | `6.3363817630495945e-15` | 6.3364 |
| `ship` | `[4/7]` | `8.8538148893157201e-15` | 8.8538 |

The continuous-extremum bounds in the two JSON files are slightly stronger
(minimum `1.3333887467657154e-15` at `2^5` and `5.3664228256360938e-15` at
`ship`). Because every admissible pair at counts 14 and 12 has a bound above
`1e-15`, no guard-admissible rational of any type at those counts meets under
the protocol. The counts 15 and 13 are therefore both attained and minimal.

## Caveat

Neither meeting fit carries a full equioscillation certificate. The targeted
Newton stopped at its iteration limit (`newton_ok` false,
`equioscillation_certificate` false) with grid/reference ratios
`1.3752977236590628` for `2^5` `[6/8]` and `1.1394203925421675` for `ship`
`[6/6]`. The count claim does not need minimax optimality of the certifying
fit: existence at the count is a constructive admissible fit scored at 60
digits with the guard and checked on 2,001, 4,001 and 20,001 nodes, and
exclusion at the count below is a de la Vallee Poussin bound, which holds for
every competitor of that type whether or not the bound-producing fit is
optimal.

The exclusion is relative to the protocol's admissible set: `m >= 1`,
`1 <= k <= 12`, and the fixed 2,001-node grid. It is not a statement about
rationals outside those limits. The independence is algorithmic: a different
solver family, different arithmetic and a different starting construction; it
shares the record's interval, target, grid and acceptance definition by
design.

## Commands

All commands were run from the scratch directory
`_clean-edge-solver-independence-2026-11-24` with `python -I`. The runtimes
are those recorded in `out/commands.txt`. Exit codes were measured by
replaying the twelve core commands in a copy of the scratch directory on
2026-11-24: all twelve exited zero and reproduced every numerical field, with
only generation timestamps and wall or round times differing.

| # | Command | Exit | Runtime |
| --- | --- | --- | --- |
| 1 | `python -I verify_scoring.py` | 0 | 23.1 s |
| 2 | `python -I remez_rational.py --test` | 0 | 1.0 s |
| 3 | `python -I sk_pairs.py --width-name "2^5" --target "rel:1e-15" --stored 15 --jobs 8 --out out/p5_s15.json` | 0 | 60.8 s |
| 4 | `python -I sk_pairs.py --width-name "2^5" --target "rel:1e-15" --stored 14 --jobs 8 --out out/p5_s14.json` | 0 | 48.2 s |
| 5 | `python -I sk_pairs.py --width-name ship --target "rel:1e-15" --stored 13 --jobs 8 --out out/ship_s13.json` | 0 | 36.7 s |
| 6 | `python -I sk_pairs.py --width-name ship --target "rel:1e-15" --stored 12 --jobs 8 --out out/ship_s12.json` | 0 | 32.9 s |
| 7 | `python -I lower_bounds.py --json out/p5_s14.json --width-name "2^5" --stored 14 --out out/lb_p5_s14.json` | 0 | 5.0 s |
| 8 | `python -I lower_bounds.py --json out/ship_s12.json --width-name ship --stored 12 --out out/lb_ship_s12.json` | 0 | 4.3 s |
| 9 | `python -I targeted_newton.py --json out/ship_s13.json --m 6 --k 6 --width-name ship --out out/targeted_ship_66.json` | 0 | ~40 s |
| 10 | `python -I targeted_newton.py --json out/p5_s15.json --m 6 --k 8 --width-name "2^5" --out out/targeted_p5_68.json` | 0 | ~35 s |
| 11 | `python -I final_checks.py` | 0 | 39.1 s |
| 12 | `python -I assemble_summary.py` | 0 | <1 s |

Exploratory runs, recorded in `out/commands.txt`, are not needed to reproduce
the result:

| Command | Exit | Runtime |
| --- | --- | --- |
| `python -I remez_rational.py --a 62.0 --width 166.72481653881133 --m 6 --k 8` | 0 | 127 s |
| `python -I diagnose_cont.py` | 0 | 527.8 s |

The first ended with the spurious `1.753192978686451e-15` local solution; its
transcript is `out/trial_p5_6_8.json`. The second is the denominator-degree
continuation; its transcript is `out/diagnose_cont.log`. Both transcripts end
with the final result line, so the recorded runs exited zero; they were not
replayed.

## Files

The seven core scripts:

| File | Contents |
| --- | --- |
| `remez_rational.py` | independent rational Remez/Newton solver, 60-digit scoring and the known-answer test |
| `sk_pairs.py` | weighted-least-squares pair sweep for one width and stored count |
| `lower_bounds.py` | de la Vallee Poussin bounds from the stored weighted-least-squares fits |
| `targeted_newton.py` | targeted Newton/Remez polish of one stored fit |
| `final_checks.py` | 60-digit scoring of the meeting fits on the 2,001-, 4,001- and 20,001-node grids |
| `assemble_summary.py` | assembles `out/summary.json` and `out/summary.txt` from the outputs |
| `verify_scoring.py` | independent re-score of the control's accepted fits from their stored coefficients |

The thirteen exploratory scripts (`diagnose_*.py`, `probe_*.py`,
`homotopy_pair.py`, `polish_pair.py`, `trial_strategies.py`) record the route
to the working initialisation. They are deposited for completeness and are
not needed to reproduce the result.

The outputs:

| File | Contents |
| --- | --- |
| `out/summary.json` | assembled evidence: inputs and digests, protocol, control cells, solver notes, counts, bounds, dense checks and conclusion |
| `out/summary.txt` | human-readable summary |
| `out/commands.txt` | the run's commands and wall times |
| `out/p5_s15.json`, `out/p5_s14.json` | `2^5` pair sweeps at stored counts 15 and 14 |
| `out/ship_s13.json`, `out/ship_s12.json` | `ship` pair sweeps at stored counts 13 and 12 |
| `out/p5_s15_68_single.json` | auxiliary single-process re-run of the `2^5` `[6/8]` fit; its `target` field reads `1e-08`, so it carries the `1e-08` record row, while the scored `e60` equals the `1e-15` value; it is not used in the summary |
| `out/lb_p5_s14.json`, `out/lb_ship_s12.json` | de la Vallee Poussin bounds for every pair at counts 14 and 12 |
| `out/targeted_p5_68.json`, `out/targeted_ship_66.json` | the two meeting fits, with 70-digit coefficients and all Newton rounds |
| `out/final_checks.json` | dense-grid scores and extrema counts |
| `out/verify_scoring.json` | control re-score on three grids |
| `out/targeted_p5_68.log`, `out/targeted_ship_66.log` | Newton round transcripts |
| `out/diagnose_cont.log`, `out/trial_strategies.log` | exploratory route transcripts |
| `out/trial_p5_6_8.json` | cold Remez transcript: a text iteration log followed by the JSON record, not valid JSON as a whole |

Four zero-byte logs in the scratch directory (`homotopy_ship_66.log`,
`polish_p5_67.log`, `polish_ship_66.log`, `trial_p5_6_8.log`) recorded no
output and are not deposited. Caches and the unused `chebfun` wheel from the
package search are also omitted.

## Identity

| Item | SHA-256 |
| --- | --- |
| `../clean-edge-2026-11-22/edge62_sweep.json` | `963facd9e18c96b547b6543c8713d5dd8294104d2776fb8d71dd30c0d40bff8b` |
| `../clean-edge-60dig-2026-11-24/cells/p5_rel_1e-15.json` | `b7dd16a6172237b250092dc6c52c0708b4514b5af790eb04bf5f821cb60ba052` |
| `../clean-edge-60dig-2026-11-24/cells/ship_rel_1e-15.json` | `70cf661f15f9f4cbc9d5665b2df40cebb65c549bf8ae7905809ade6cc845455b` |
| `out/summary.json` | `218e63e2de8af56a118e58ebbce0dcbd65b77d949538b2d6de65ab04d8ec9f0e` |
| `out/summary.txt` | `17b85c5b09a17b5d765976fa8b2aae8182036c5945af4f4eb2619f79a3ef1b90` |
| `out/commands.txt` | `17d1af76a3f156043b94682ab9b22005fa43535ff2b329e8c964a36c9b1cd6a6` |
| `out/p5_s15.json` | `af73598cd24926ccbe6772b984d02d449d689dea1570990533e48fd00f17872b` |
| `out/p5_s14.json` | `9b8d98ed75b7e7fac839857d5d5fa690c295dc831748d79611b834fb4584d7cf` |
| `out/ship_s13.json` | `34a773c2e22ac70883e3b0a55d4587a953ad75024de816f99c6d8182b5db8661` |
| `out/ship_s12.json` | `d9652cd751e8f2ee7ca08e075e898b964e29f7320dbf297c7ca6c7e37448bc40` |
| `out/lb_p5_s14.json` | `26353b7a69f0a97194d36fb5061b17acc4393564af36dc8f8b7052bfd34eba39` |
| `out/lb_ship_s12.json` | `f9c2e4a8ff1a7eb3b33963d9cbdb0af6b5f095a7a52bbbb40dfd3681208efbd8` |
| `out/targeted_p5_68.json` | `4dbe01bac07f6808cad239bf1b4f24408994288bba6160c941efd0c5d3f90fe0` |
| `out/targeted_ship_66.json` | `754bf370ce3b31f13317433196ad7efc73bd09a7955fcc487a7e816b4fb180af` |
| `out/final_checks.json` | `1d19045b69c977a426646d67192739faf69fa045dbd1ecac2164b07a1a3be965` |
| `out/verify_scoring.json` | `0d0a39228828cbb58ad38512a1c13ecd6d9b6c87a96fbcdcf943f900f504372f` |
| `remez_rational.py` | `525e41282508fc0b64cd2d4e9d1ae9e90b399bccf1a970eea7dfc621ddc0fcd3` |
| `sk_pairs.py` | `afa64696596279f22bf32d35361682127a88cb3fcff1a87787fe9764fca0f306` |
| `lower_bounds.py` | `3ebb89715ba4e984de33f57b0b66bdb0d6e9c4a38d335a4692766ffe09d6b6b3` |
| `targeted_newton.py` | `f3177dea37e52d98879cadbbf8ed8b0bf7b04491042666cd5b6d4b6b6e674c98` |
| `final_checks.py` | `eb5ce97c0cb6e7fc2b1b241bfa8325f0a04e2a81c12fb1848e8bbf129c24b6da` |
| `assemble_summary.py` | `183d54662f4751c615c255600f8a82bb6cb1d80fec96f4ab49dd8090d2c46efd` |
| `verify_scoring.py` | `b03f65410a0944f3862fde7f8e140c9a7de9c69b169a865f3a4ecf5588ff94e4` |
| `environment.txt` | `e45273bdce480505d9bfac643f95ca0085778a78d518a6ec1f226963841bf231` |

The first three rows are the binary64 clean-edge record and the two 60-digit
control cells, all deposited in this repository. The scripts' digests are
those of the deposited copies.

## Environment

`environment.txt` records the environment captured from the scratch directory:
Windows 11 (`10.0.26200`), Python 3.14.7, `mpmath` 1.4.1 with its `gmpy`
backend (`gmpy2` 2.3.1), `numpy` 2.5.2, `scipy` 1.18.1 and `baryrat` 2.1.2.
The independent arm imports only `mpmath` beside the standard library;
`numpy`, `scipy` and `baryrat` belong to the surrounding control environment
and are not imported by any deposited script.

## Reproducing

The scripts carry the measurement machine's absolute input paths. Point them
at the deposit's repository siblings before a rerun:

- `REC` in `sk_pairs.py`, `lower_bounds.py`, `targeted_newton.py`,
  `final_checks.py`, `assemble_summary.py`, `homotopy_pair.py` and
  `polish_pair.py` at `../clean-edge-2026-11-22/edge62_sweep.json`;
- `CELLS` in `verify_scoring.py`, `assemble_summary.py`, `diagnose_init.py`
  and `diagnose_p5.py` at `../clean-edge-60dig-2026-11-24/cells`.

The deposited record and cells match the measurement inputs byte for byte
(digests above), so a rerun reads the same data. Python 3.10 or newer with
`mpmath` is enough. Run from this directory so the `out/...` arguments land
in `out/`:

```sh
python -I remez_rational.py --test
python -I verify_scoring.py
python -I sk_pairs.py --width-name "2^5" --target "rel:1e-15" --stored 15 --jobs 8 --out out/p5_s15.json
python -I sk_pairs.py --width-name "2^5" --target "rel:1e-15" --stored 14 --jobs 8 --out out/p5_s14.json
python -I sk_pairs.py --width-name ship --target "rel:1e-15" --stored 13 --jobs 8 --out out/ship_s13.json
python -I sk_pairs.py --width-name ship --target "rel:1e-15" --stored 12 --jobs 8 --out out/ship_s12.json
python -I lower_bounds.py --json out/p5_s14.json --width-name "2^5" --stored 14 --out out/lb_p5_s14.json
python -I lower_bounds.py --json out/ship_s12.json --width-name ship --stored 12 --out out/lb_ship_s12.json
python -I targeted_newton.py --json out/ship_s13.json --m 6 --k 6 --width-name ship --out out/targeted_ship_66.json
python -I targeted_newton.py --json out/p5_s15.json --m 6 --k 8 --width-name "2^5" --out out/targeted_p5_68.json
python -I final_checks.py
python -I assemble_summary.py
```

A fresh run reproduces the numerical fields. `summary.json`, `summary.txt`
and the pair sweeps differ in their generation timestamps and recorded wall
times; the identity table applies to the deposited files.
