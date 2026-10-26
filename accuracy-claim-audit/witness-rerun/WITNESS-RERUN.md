# Witness rerun for the two false-sentence findings

Bounded re-derivation of the two witnesses reported in the accuracy audit of
Section `sec:audit`. The rerun was made on 26 October 2026 from public,
pinned sources, on one Windows 11 machine with WSL2 for the Linux build. Both
documented sentences are re-tested against an mpmath reference at 60
significant digits. The witness directory is deposited as
`accuracy-claim-audit/witness-rerun/`; the deposit is the bounded set
described in Section 4, with the upstream sources pinned by hash rather than
redistributed. The commands, logs, every argument evaluated, and every output
are there. This report carries the tables and the verdicts. No other library
was substituted for either claim.

The two verdicts in brief. The Boost.Math Airy sentence fails at both
arguments the audit names: `airy_ai` at `x = 35.9` is 121.395340 eps from the
reference, and `airy_bi` at the audit's argument behind the printed
`x = 78.1252` is 414.110896 eps, reproducing the audit's 121.395 and 414.111
to the printed three decimals. The literal four-decimal double 78.1252 also
breaches the sentence, at 149.324772 eps, so the verdict does not depend on
reconstructing the audit grid. crlibm's `sinh` entry also fails as reported:
at `0x408633ce8fb9f87e`, `sinh_rd` and `sinh_rz` return `+inf` where IEEE 754
clause 7.4 requires `+DBL_MAX`, and at the negated input `sinh_ru` and
`sinh_rz` return `-inf` where `-DBL_MAX` is required. `sinh_rn` is correct at
that input because the exact value lies above the round-to-nearest overflow
midpoint; the paper's wording does not claim otherwise.

## 1. Sources, tools, and commands

### 1.1 Pinned sources

| Source | Pin | SHA-256 | Version stamp |
| --- | --- | --- | --- |
| Boost 1.92.0 | `https://archives.boost.io/release/1.92.0/source/boost_1_92_0.tar.gz` | `c4a3b310ddd2472416e091067166b0713be97c63f38c212c484ada022fd296ce` | `BOOST_VERSION 109200`, `BOOST_LIB_VERSION "1_92"` |
| crlibm | `https://github.com/taschini/crlibm/archive/eb3063791aa75bc9705b49283bf14250465220a7.tar.gz` | `79d626f9515fccdbc724e8b74cbf848331317b1891ef22aaaf6a195e20bbbfe5` | commit `eb3063791aa75bc9705b49283bf14250465220a7`, `VERSION` 1.0beta5 |

The Boost hash matches the release site's own `.json` descriptor for 1.92.0.
The crlibm archive is the GitHub mirror of the INRIA Forge release, at the
commit the audit used. The archives and extracted trees are not redistributed;
the pins fix their bytes. The documented sentences were quoted from the
pinned sources at run time; the excerpts and their locations are in
`raw/doc_excerpts.txt`.

### 1.2 Toolchain and flags

| Item | Version or setting |
| --- | --- |
| Windows compiler | Microsoft C/C++ Optimizing Compiler 19.51.36260 for x64 (Visual Studio 18 Community) |
| Boost witness flags | `/nologo /std:c++17 /EHsc /O2 /W3 /MD /fp:precise` |
| Boost mode | header-only, `boost/math/special_functions/airy.hpp` |
| Python | 3.14.7 (MSC v.1944, 64-bit) |
| mpmath | 1.4.1, `mp.dps = 60` (80 for the threshold context) |
| Linux compiler | gcc 15.2.0 (Ubuntu 15.2.0-16ubuntu1), x86-64 |
| crlibm flags | `-std=c99 -Wall -msse2 -mfpmath=sse -O2` |
| crlibm configure | `CFLAGS="-O2" ./configure --enable-sse2 --disable-mpfr --disable-double-extended` |
| Driver header at run | `FLT_EVAL_METHOD=0`, `sizeof(double)=8`, `fegetround()=0` |

The Python scripts were run with `python -I`, which hides user
`site-packages`; mpmath is installed in the interpreter's main
`site-packages` here.

The audit's README describes crlibm as built in the non-double-extended
configuration but its listed configure command omits
`--disable-double-extended`. `configure.ac` gives that branch as the default,
and on this host the default selected `log-de.c`. The rerun passes the flag
and selects `log.c`. `csh_fast.c`, the file that implements `sinh`, has no
double-extended conditional and is compiled from the same source either way;
the flag changes the `log` implementation, not the `sinh` code under test.

### 1.3 Commands and exit codes

From the root of the witness directory. The full record, including the log
redirections and the two second runs, is in `logs/commands.txt`.

| Command | Exit |
| --- | --- |
| `python -I scripts/make_boost_inputs.py` | 0 |
| `cmd /c scripts\build_boost.bat > logs\boost_build.log 2>&1` | 0 |
| `raw/boost_airy_witness.exe raw/boost_airy_inputs.txt > raw/boost_airy_raw.csv` | 0 |
| `python -I scripts/reference_boost.py > logs/reference_boost.log 2>&1` | 0 |
| `python -I scripts/make_crlibm_inputs.py` | 0 |
| `bash scripts/build_crlibm.sh > logs/crlibm_build.log 2>&1` (in WSL, final attempt) | 0 |
| `raw/crlibm_sinh_witness raw/crlibm_sinh_inputs.txt > raw/crlibm_sinh_raw.csv` (in WSL) | 0 |
| `python -I scripts/reference_crlibm.py > logs/reference_crlibm.log 2>&1` | 0 |
| `python -I scripts/crlibm_threshold_context.py` | 0 |
| `cmp raw/boost_airy_raw.csv logs/boost_airy_raw_rerun.csv` | 0 |
| `cmp raw/crlibm_sinh_raw.csv logs/crlibm_sinh_raw_rerun.csv` | 0 |

The first crlibm build attempt exited 1 at the driver compile because the
driver used `FLT_EVAL_METHOD` without including `<float.h>`; the library
itself had built. The header was added and the build rerun. The same attempt
had also selected the double-extended branch, which
`--disable-double-extended` now prevents. The deposited build log is the
successful run, and it shows `log.c` and the flags above. The two `cmp`
commands produced no output; the deposited `logs/*_raw_rerun.csv` files are
byte-identical to the first runs.

## 2. Claim 1: Boost.Math's Airy sentence

### 2.1 The documented sentence

The sentence appears in the Boost 1.92.0 documentation on the `airy_ai` and
`airy_bi` pages. The source is `libs/math/doc/sf/airy.qbk`, lines 52-53 (Ai)
and 106-107 (Bi), and again on the two prime pages:

```text
In general though, the relative error is low (less than 100 [epsilon]) for /x > 0/ while
only the absolute error is low for /x < 0/ as the following error plot illustrates:
```

The rendered pages print `[epsilon]` as the Greek epsilon. The paper's
quotation, "less than 100 eps", is this sentence; epsilon is the binary64
format epsilon, 2^-52. The domain is the stated qualifier `x > 0`, with no
upper bound.

### 2.2 Method

The witness driver reads 64-bit patterns, evaluates
`boost::math::airy_ai(x)` and `boost::math::airy_bi(x)` in `double`, and
prints the result bits together with a decimal rendering. The bits are
authoritative; the decimal is only for reading. It was built with MSVC
Release at `/O2` and `/fp:precise`, as the audit used MSVC. The environment
header in the raw CSV records `BOOST_VERSION 109200`, `_MSC_VER 1951`,
`_M_FP_PRECISE 1`, `sizeof(double) 8`, and `FLT_EVAL_METHOD 0`.

The reference is mpmath 1.4.1 at 60 significant decimal digits, evaluated at
the exact binary64 value of each input, never at its printed decimal. The
error is reported in units of eps = 2^-52:

```text
eps_units = |library - reference| / |reference| / 2^-52
```

The input list has 2349 distinct arguments: the named arguments and their
immediate neighbours, dense windows of step 0.001 around each named point,
coarse windows, and the two audit grids reconstructed from the row table
(`linear:65` over `[0.001, 100]` and `linear:33` over `[30, 35.9]`).

Three checks support the reference. At `x = 35.9` the Airy ODE residual is
`0.0` for Ai and `8.15161e-62` for Bi relative; the Wronskian `Ai Bi' - Ai' Bi`
is `0.318309886183790671537767526745`, which is `1/pi` to a relative residual
of `1.22189e-61`; and every named point and window maximum rounds to the same
binary64 at 60 and at 100 digits. The quoted output is in
`raw/boost_independent_check.txt`.

### 2.3 Measured errors at the named arguments

| Function | Argument | Input bits | Error (eps) | Audit/paper |
| --- | --- | --- | --- | --- |
| `airy_ai` | 35.9 | `4041f33333333333` | 121.395340 | 121.395 |
| `airy_ai` | 78.1252 (literal double) | `4053880346dc5d64` | 148.297779 | — |
| `airy_ai` | 78.12521875 (audit grid i=50) | `4053880395810625` | 414.878957 | — |
| `airy_ai` | 78.13 | `40538851eb851eb8` | 247.193968 | — |
| `airy_bi` | 35.9 | `4041f33333333333` | 120.848885 | — |
| `airy_bi` | 78.1252 (literal double) | `4053880346dc5d64` | 149.324772 | — |
| `airy_bi` | 78.12521875 (audit grid i=50) | `4053880395810625` | 414.110896 | 414.111 |
| `airy_bi` | 78.13 | `40538851eb851eb8` | 244.728194 | — |
| `airy_bi` | 88.258434442270058 (gate context) | `4056108a309d6d6f` | 529.766809 | 529.767 |

The exact-value comparison at the load-bearing points:

| Function | x (bits) | Library value | mpmath reference (60 digits) | Relative error | eps |
| --- | --- | --- | --- | --- | --- |
| `airy_ai` | 35.9 (`4041f33333333333`) | `6.072983813676663e-64` | `6.07298381367682676012450295236e-64` | `2.6955180215907676e-14` | 121.395340 |
| `airy_bi` | 78.1252 (`4053880346dc5d64`) | `1.6185318944921704e+199` | `1.61853189449222407526347701457e+199` | `3.315675998556602e-14` | 149.324772 |
| `airy_bi` | 78.12521875 (`4053880395810625`) | `1.618800056163607e+199` | `1.61880005616375578780131223317e+199` | `9.1951090213725588e-14` | 414.110896 |
| `airy_ai` | 78.12521875 (`4053880395810625`) | `1.1123243383936614e-201` | `1.11232433839355888728847712376e-201` | `9.2121634155000135e-14` | 414.878957 |

Two dense-window maxima locate the same failure away from the named points:

```text
airy_ai  window [35.4, 36.4]: max 127.149510 eps at x=36.282 (404224189374bc6a)
airy_bi  window [77.6, 78.6]: max 441.685988 eps at x=77.764 (405370e560418937)
```

The `bi` maximum in the window around the printed `78.1252` is 441.686 eps,
at 77.764. The failure is not a knife-edge at one double.

### 2.4 The argument behind the printed "78.1252"

The audit's row table records the row `A2` as `linear:65` over
`[0.001, 100]`. A 65-point linear grid has spacing `99.999/64`, and point
`i = 50` is

```text
0.001 + 50 * (99.999 / 64) = 78.12521875
```

whose binary64 pattern is `4053880395810625`. At that double the rerun
measures `airy_bi = 414.110896 eps`, matching the audit's `414.111` to the
printed three decimals. The paper prints the argument at four decimals,
`78.1252`, which is the same grid point rounded for display; the literal
four-decimal double is a different binary64 value, `4053880346dc5d64`, and
there the measured error is `149.324772 eps`.

The difference is not an error in either measurement. The Airy error is a
sawtooth in the last bits of the argument: the implementation forms
`(2/3)x^{3/2}` in binary64, and a one-ulp change in `x` moves the exponent of
the exponential by a comparable amount, so the relative error changes by
tens of eps between adjacent doubles. The audit evaluated the grid double
and displayed it at four decimals; the rerun evaluates both the grid double
and the literal label. Both are in the deposited input list, and both
breach the sentence.

The paper's later sentence, "At the widest sampled argument, x = 78.13",
refers to a different double and a different part of the argument: the rerun
measures 247.193968 eps for `airy_ai` and 244.728194 eps for `airy_bi` at
78.13, consistent with the paper's half-ulp floor discussion there. The
414.111 audit cell is the A2 grid point `78.12521875`, not 78.13; the two
arguments are 0.0048 apart and their errors differ because the sawtooth has
already moved.

### 2.5 Raw output

The driver's environment header and the three relevant rows, quoted
verbatim from `raw/boost_airy_raw.csv`:

```text
# BOOST_VERSION=109200 BOOST_LIB_VERSION=1_92
# _MSC_VER=1951 _MSC_FULL_VER=195136260
# _M_FP_PRECISE=1
# sizeof(double)=8 FLT_EVAL_METHOD=0
input_bits,x,airy_ai_bits,airy_ai,airy_bi_bits,airy_bi
4041f33333333333,35.899999999999999,32cffa5db1483540,6.0729838136766631e-64,4cbb38144b7056e4,4.3739342516351482e+61
4053880346dc5d64,78.125200000000007,1635ccda11aa3a7b,1.1125087640638141e-201,694b10bcab441d94,1.6185318944921704e+199
4053880395810625,78.125218750000002,1635cbed3a28cf03,1.1123243383936614e-201,694b11e28c737f82,1.6188000561636069e+199
```

The corresponding reference rows, read from
`raw/boost_airy_reference.csv` and laid out with the columns (bits, function,
result, reference, eps):

```text
4041f33333333333  airy_ai  6.072983813676663e-64     6.07298381367682676012450295236e-64   121.395339576
4041f33333333333  airy_bi  4.373934251635148e+61     4.37393425163503084973065827892e+61   120.848885103
4053880346dc5d64  airy_bi  1.6185318944921704e+199   1.61853189449222407526347701457e+199   149.324771916
4053880395810625  airy_bi  1.618800056163607e+199    1.61880005616375578780131223317e+199   414.110895623
```

### 2.6 Verdict

The documented sentence fails, as the paper reports, and the rerun reproduces
both witnesses to the printed precision. At `airy_ai`, `x = 35.9`, the error
is 121.395340 eps against the stated "less than 100 eps", a 21.4 per cent
breach of the figure; the audit's 121.395 reproduces exactly at the literal
argument. At `airy_bi`, the audit's argument is the grid double
`78.12521875`, displayed in the paper as `78.1252`; there the error is
414.110896 eps, reproducing the audit's 414.111, a 314.1 per cent breach. The
literal four-decimal double `78.1252` is a different value and gives
149.324772 eps, a 49.3 per cent breach, so the sentence is false at the
argument as printed even if the grid reconstruction is set aside. The
nearest dense-window maximum, 441.686 eps at 77.764, is higher still. The
sentence's domain is `x > 0` without an upper bound, and these arguments lie
inside it; the claim that the relative error is low, under the natural
reading of "less than 100 eps", is therefore false at both audited functions.
The hedge "in general" leaves the sentence unfalsifiable only if it is read
as having no bound at all; the paper already states that reading it as a
bound is an interpretation, and under that reading the figure is breached by
the margins above. The rerun measures mpmath at 60 digits with ODE,
Wronskian and 60-vs-100 digit checks; it does not re-run the
Boost.Multiprecision arm, so the paper's separate 2.7 eps cross-reference
agreement is not re-measured here.

## 3. Claim 2: crlibm's sinh entry

### 3.1 The documented sentence

The crlibm 1.0beta5 manual, `docs/latex/0_getting-started.tex`, line 149,
lists in Table 1:

```text
sinh & sinh\_rn & sinh\_ru & sinh\_rd & sinh\_rz & complete & complete (paper)\\ \hline
```

The table's own definition of the first state column, lines 173-177, is:

```text
\item The first indicates the state of the search for worst cases
  for correct rounding \cite{LMT98,Lef2000}. If it indicates
  ``complete'', it means that the function is guaranteed to return
  correct rounding on its whole floating-point input range.
```

One counterexample to correct rounding over that range refutes the entry.

### 3.2 Method

crlibm 1.0beta5 was built from the pinned mirror in WSL with gcc 15.2.0 in
the non-double-extended configuration, with `-O2` and the SSE2 flags shown
in Section 1.2. The driver calls `sinh_rn`, `sinh_ru`, `sinh_rd` and
`sinh_rz` on each input and prints the four result bit patterns; bits are
authoritative. It was compiled with gcc at `-O2 -std=c99 -msse2 -mfpmath=sse`
and linked against the built `libcrlibm.a`.

The input list has 422 distinct arguments: the named double
`0x408633ce8fb9f87e`, 64 neighbours on each side, the negative mirror
`0xc08633ce8fb9f87e` and its neighbours, a scan of the positive and negative
overflow tail in the form the audit's `m3_defects.py` uses, `DBL_MAX`, two
points below the boundary, and 1.0 as a sanity check.

The reference is mpmath 1.4.1 at 60 significant digits. The exact value is
kept as an exact mpf, converted to a `Fraction` when its magnitude is at most
`DBL_MAX`, and rounded to binary64 under each direction with the clause 7.4
overflow rule:
a magnitude above `DBL_MAX` overflows; under round-to-nearest it carries to
infinity only at or above the midpoint `T_RN = 2^1024 - 2^970`, while
toward-zero and toward-minus-infinity carry a positive overflow to `DBL_MAX`
and toward-plus-infinity carries it to `+inf`. The expected bits for all four
modes are in `raw/crlibm_sinh_reference.csv`.

The reference's own checks: the ODE residual `sinh'' - sinh` is `0.0`
relative, `cosh(1/3)^2 - sinh(1/3)^2 - 1` is `0.0`, `log(sinh(x))` equals
709.782712893384031803813853119, `asinh(exact sinh) - x` is `0.0`, and the
rounding decision is the same at 60 and at 100 digits for every input and
mode (0 differences). The quoted output is in
`raw/crlibm_sinh_independent_check.txt`.

### 3.3 The exact value and the overflow thresholds

At `0x408633ce8fb9f87e` the exact double is
710.4758600739439771132310852408409118652 and the exact hyperbolic sine,
computed at 80 digits, is

```text
sinh(x) = 1.7976931348623787561026786514116073629845935e+308
```

The two thresholds are

```text
DBL_MAX = 2^1024 - 2^971 = 1.797693134862315708145274237317043567981e+308
T_RN    = 2^1024 - 2^970 = 1.797693134862315807937289714053034150799e+308
```

The exact value is above both. Relative to `DBL_MAX` it is larger by
`3.50715904632e-14`; relative to `T_RN` it is larger by `3.5016079312e-14`.
Equivalently, `asinh(DBL_MAX) = 710.475860073943942041640622032...` and
`asinh(T_RN) = 710.475860073943942097151773263...`, both below the input
double. The input is therefore the first double at which `sinh` overflows
under round-to-nearest, and also the first at which it overflows under every
directed mode. This is why `sinh_rn` returning `+inf` is correct and why
`sinh_rd` and `sinh_rz` must return `+DBL_MAX`, not `+inf`. The numbers are in
`raw/crlibm_threshold_context.txt`.

### 3.4 Measured results at the witness and its neighbours

| Input bits | Mode | Result bits | Correct? | Required |
| --- | --- | --- | --- | --- |
| `408633ce8fb9f87d` (below) | `sinh_rd` | `7feffffffffffd3a` | yes | same |
| `408633ce8fb9f87e` (named) | `sinh_rn` | `7ff0000000000000` | yes | `+inf` (above `T_RN`) |
| `408633ce8fb9f87e` (named) | `sinh_ru` | `7ff0000000000000` | yes | `+inf` |
| `408633ce8fb9f87e` (named) | `sinh_rd` | `7ff0000000000000` | **no** | `7fefffffffffffff` (`+DBL_MAX`) |
| `408633ce8fb9f87e` (named) | `sinh_rz` | `7ff0000000000000` | **no** | `7fefffffffffffff` (`+DBL_MAX`) |
| `408633ce8fb9f87f` (above) | `sinh_rd` | `7fefffffffffffff` | yes | same |
| `c08633ce8fb9f87e` (mirror) | `sinh_rn` | `fff0000000000000` | yes | `-inf` |
| `c08633ce8fb9f87e` (mirror) | `sinh_ru` | `fff0000000000000` | **no** | `ffefffffffffffff` (`-DBL_MAX`) |
| `c08633ce8fb9f87e` (mirror) | `sinh_rd` | `fff0000000000000` | yes | `-inf` |
| `c08633ce8fb9f87e` (mirror) | `sinh_rz` | `fff0000000000000` | **no** | `ffefffffffffffff` (`-DBL_MAX`) |

The exact comparisons for the same points:

| Input bits | Exact `sinh` | `rn` expected / got | `ru` expected / got | `rd` expected / got | `rz` expected / got |
| --- | --- | --- | --- | --- | --- |
| `408633ce8fb9f87d` | `1.79769313486217438205498230058e+308` | `7feffffffffffd3b` / same | `7feffffffffffd3b` / same | `7feffffffffffd3a` / same | `7feffffffffffd3a` / same |
| `408633ce8fb9f87e` | `1.79769313486237875610267865141e+308` | `7ff0000000000000` / same | `7ff0000000000000` / same | `7fefffffffffffff` / `7ff0000000000000` | `7fefffffffffffff` / `7ff0000000000000` |
| `408633ce8fb9f87f` | `1.79769313486258313015037502548e+308` | `7ff0000000000000` / same | `7ff0000000000000` / same | `7fefffffffffffff` / same | `7fefffffffffffff` / same |
| `c08633ce8fb9f87e` | `-1.79769313486237875610267865141e+308` | `fff0000000000000` / same | `ffefffffffffffff` / `fff0000000000000` | `fff0000000000000` / same | `ffefffffffffffff` / `fff0000000000000` |

Over the whole sampled set, the mismatch counts are:

| Mode | Mismatches | Of |
| --- | --- | --- |
| `sinh_rn` | 0 | 422 |
| `sinh_ru` | 1 | 422 |
| `sinh_rd` | 1 | 422 |
| `sinh_rz` | 2 | 422 |

The only mismatches are at two inputs: the positive named double, where
`sinh_rd` and `sinh_rz` fail, and the negative mirror, where `sinh_ru` and
`sinh_rz` fail. The sampled inputs
where a routine returned `+inf` are 146 sampled positive inputs spanning the
named double up to `DBL_MAX` for `sinh_rn` and `sinh_ru`, none wrong; exactly
one input for `sinh_rd`, the named double, wrong; and one positive input for
`sinh_rz`, also the named double, wrong. By the arithmetic and the source's
fast-path guard, every positive double above the named one returns `+inf`
from `sinh_rn` and `sinh_ru`, and that is correct; `sinh_rd` and `sinh_rz`
return `DBL_MAX` from the next double upward, also correct. The full list of
sampled inputs and exact values is in `raw/crlibm_sinh_reference.csv`; at
`DBL_MAX` the exact `sinh` is far above `T_RN`, and all four modes return
their clause 7.4 results.

### 3.5 The mechanism

The four routines share a guard whose comparison is strict:

```c
/* csh_fast.c, line 671 (sinh_rd); the same form at line 619 (sinh_ru) */
if (y.d > max_input_csh.d) { /* out of range */
```

At the double exactly equal to `max_input_csh.d` the guard is not taken, so
`sinh_rd` enters the accurate path. That path rounds the significand
downwards, but its final step restores the exponent with a binary64
multiplication:

```c
/* csh_fast.c, line 706 (sinh_rd); the same form at line 656 (sinh_ru) */
res.i[HI] += (exponent-11) << 20;
return 1024. * res.d;
```

The multiplication runs in the caller's rounding mode, which is
round-to-nearest here, so for a true value above `T_RN` it overflows to
`+inf` and discards the directed rounding just performed. `sinh_rz` is a
dispatcher: it calls `sinh_rd` for positive inputs and `sinh_ru` for negative
ones, so it inherits both faults. At the next double upward the strict guard
does catch the input and returns `DBL_MAX` directly, which is why the defect
is exactly one double per sign; at the next double downward the true value is
below `DBL_MAX` and no overflow occurs. The code comment on the exponent
adjustment says the scaling "will cater for overflow", but the multiplication
does so in round-to-nearest, which is the wrong direction for `rd` and, on
the negative side, for `ru`.

### 3.6 Raw output

The driver header and the named input, mirror, and neighbours, quoted
verbatim from `raw/crlibm_sinh_raw.csv`:

```text
# FLT_EVAL_METHOD=0 sizeof(double)=8 sizeof(long double)=16 fegetround_at_start=0
input_bits,x,sinh_rn_bits,sinh_rn,sinh_ru_bits,sinh_ru,sinh_rd_bits,sinh_rd,sinh_rz_bits,sinh_rz
408633ce8fb9f87d,710.47586007394386,7feffffffffffd3b,1.7976931348621744e+308,7feffffffffffd3b,1.7976931348621744e+308,7feffffffffffd3a,1.7976931348621742e+308,7feffffffffffd3a,1.7976931348621742e+308
408633ce8fb9f87e,710.47586007394398,7ff0000000000000,inf,7ff0000000000000,inf,7ff0000000000000,inf,7ff0000000000000,inf
408633ce8fb9f87f,710.47586007394409,7ff0000000000000,inf,7ff0000000000000,inf,7fefffffffffffff,1.7976931348623157e+308,7fefffffffffffff,1.7976931348623157e+308
c08633ce8fb9f87e,-710.47586007394398,fff0000000000000,-inf,fff0000000000000,-inf,fff0000000000000,-inf,fff0000000000000,-inf
```

The corresponding reference rows, read from
`raw/crlibm_sinh_reference.csv` and laid out one field per line:

```text
408633ce8fb9f87e  710.475860073944  exact 1.79769313486237875610267865141e+308
  rn expected 7ff0000000000000  got 7ff0000000000000  match True
  ru expected 7ff0000000000000  got 7ff0000000000000  match True
  rd expected 7fefffffffffffff  got 7ff0000000000000  match False
  rz expected 7fefffffffffffff  got 7ff0000000000000  match False
c08633ce8fb9f87e -710.475860073944  exact -1.79769313486237875610267865141e+308
  ru expected ffefffffffffffff  got fff0000000000000  match False
  rz expected ffefffffffffffff  got fff0000000000000  match False
```

### 3.7 Verdict

The documented "complete" entry fails, exactly as the paper reports, at the
argument the paper names. At `0x408633ce8fb9f87e`, `sinh_rd` and `sinh_rz`
return `+inf` where IEEE 754 clause 7.4 requires `+DBL_MAX`
(`7fefffffffffffff`); at the negated input `0xc08633ce8fb9f87e`, `sinh_ru`
and `sinh_rz` return `-inf` where `-DBL_MAX` (`ffefffffffffffff`) is
required. The exact value at the named input,
`1.7976931348623787561e308`, is above `DBL_MAX` by `3.51e-14` relative and
above the round-to-nearest midpoint `T_RN` by `3.50e-14` relative, so
`sinh_rn`'s `+inf` is correct; the fault is the directed modes, and the
paper's wording is precise about which routines fail. The defect is exactly
one double per sign: the other 420 sampled arguments, including the immediate
neighbours on both sides, match in all four modes. The mechanism is the strict `>` guard at
`csh_fast.c:671` (and `:619` for `sinh_ru`), which lets the boundary double
reach an accurate path whose final binary64 multiplication overflows under
round-to-nearest, undoing the directed rounding. The "complete" entry claims
correct rounding over the whole floating-point input range; one double per
sign is enough to refute it. The rerun uses mpmath at 60 digits with exact
rational rounding and the clause 7.4 rule; it does not re-run the audit's
MPFR or glibc routes, so those two independent confirmations are not
re-measured here.

## 4. Reproduction, determinism, and the deposit

The witness directory is deposited in full apart from the two pinned
upstream sources and the rebuildable binaries. `README.md` in it is a
stranger-followable reproduction: it gives the download URLs and SHA-256
values, the extraction commands, the exact build flags, the input generation,
the run commands, the reference commands, and the expected headline numbers.
`logs/commands.txt` records the command lines and observed exit codes;
`logs/environment.txt` records the tool versions, flags and hashes.
`raw/doc_excerpts.txt` quotes the two documented sentences from the pinned
sources.

Both witness programs were run a second time and compared byte for byte with
`cmp`. The comparisons exited 0 with no output; the second runs are deposited
as `logs/boost_airy_raw_rerun.csv` and `logs/crlibm_sinh_raw_rerun.csv`. The
reference scripts are deterministic given the raw CSV.

The bounded deposit excludes the two upstream archives, the extracted Boost
and crlibm trees, and the built witness binaries and object file. Those are
third-party sources, they are large (about 199 MB extracted and 237 MB of
archives), and the paper states that crlibm's source is not redistributed.
Their bytes are fixed by the SHA-256 pins in Section 1.1, and rebuilding is
the first two commands of each reproduction section. The copy of the witness
directory at `accuracy-claim-audit/witness-rerun/` is this bounded set, with
a `MANIFEST.sha256` listing every included file.

No git command was run during the rerun, and no existing file in the paper
tree was modified. This report is the only new file in the paper tree.

## 5. Limits of the rerun

The sample is bounded: 2349 Boost arguments and 422 crlibm arguments. It is
large enough to reproduce the audit's witnesses and to show they are not
isolated rounding accidents, but it is not an exhaustive scan of either
domain, and it does not claim to find the worst case over all doubles. The
Boost grid reconstruction is the one inferred step: the row table gives
`linear:65` and the interval, and the rerun reconstructs point `i = 50` as
`78.12521875`, which reproduces the audit's 414.111. The literal four-decimal
argument is also deposited and also breaches the sentence, so the verdict
stands without the reconstruction; only the exact reproduction of the
audit's figure depends on it. The Boost rerun measures mpmath at 60 digits
and does not re-run the Boost.Multiprecision arm, so the paper's separate
2.7 eps cross-reference agreement is outside this rerun. The crlibm rerun
measures mpmath at 60 digits with exact rational rounding and does not
re-run the audit's MPFR or glibc routes. The crlibm build adds
`--disable-double-extended` to match the audit's described configuration;
the audit's listed command omits it, and `csh_fast.c` has no conditional on
that branch, so the `sinh` code under test is the same source.
