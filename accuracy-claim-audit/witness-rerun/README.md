# Witness rerun: two documented accuracy sentences

This directory re-derives, in a bounded and fully deposited run, the two
false-sentence witnesses reported in the Paper B accuracy audit:

1. **Boost.Math 1.92.0, Airy functions.** The sentence on the `airy_ai` and
   `airy_bi` documentation pages says the relative error is "less than 100 ε"
   for `x > 0`. The audit reports breaches at `airy_ai`, `x = 35.9`, 121.395 ε,
   and at `airy_bi`, `x = 78.1252`, 414.111 ε.
2. **crlibm 1.0beta5, `sinh`.** Table 1 of the manual marks the `sinh` entry
   "complete", which its own definition explains as correct rounding over the
   whole floating-point input range. The audit reports `sinh_rd` and
   `sinh_rz` returning `+inf` at `0x408633ce8fb9f87e`, where IEEE 754
   clause 7.4 requires `+DBL_MAX`, with the mirrored failure in `sinh_ru` and
   `sinh_rz` at the negated input.

Everything needed to reproduce both measurements is here except the two
upstream source archives, which are pinned by URL and SHA-256 below. The
verdicts and the tables of measured values are in `WITNESS-RERUN.md` in this
directory; the paper cites the same report.

## Contents

| Path | What it is |
| --- | --- |
| `WITNESS-RERUN.md` | the report: sources and tools, the commands and exit codes, both verdicts, the exact values and thresholds, and the limits of the rerun |
| `src/boost_airy_witness.cpp` | MSVC driver: reads 64-bit patterns, prints `airy_ai` and `airy_bi` bits and decimals |
| `src/crlibm_sinh_witness.c` | gcc driver: reads 64-bit patterns, prints `sinh_rn`, `sinh_ru`, `sinh_rd`, `sinh_rz` bits and decimals |
| `scripts/make_boost_inputs.py` | writes the Boost argument list and labels |
| `scripts/reference_boost.py` | mpmath 60-digit reference, eps errors, and independent checks |
| `scripts/build_boost.bat` | MSVC Release build, `/O2 /fp:precise`, header-only Boost |
| `scripts/make_crlibm_inputs.py` | writes the crlibm argument list: the witness, its neighbours, the mirror, and a coarse scan of the overflow tail |
| `scripts/build_crlibm.sh` | autotools build of crlibm in the non-double-extended configuration, then the driver |
| `scripts/reference_crlibm.py` | mpmath 60-digit reference with the clause 7.4 overflow rule, and independent checks |
| `scripts/crlibm_threshold_context.py` | mpmath 80-digit comparison of the exact `sinh` with `DBL_MAX` and the round-to-nearest midpoint |
| `raw/*inputs*` | every argument evaluated, with labels |
| `raw/*raw.csv` | the drivers' raw output, one row per argument |
| `raw/*reference.csv` | the reference value, the rounded result expected in each mode, and the measured error |
| `raw/*independent_check.txt` | reference checks quoted in the report |
| `raw/crlibm_threshold_context.txt` | the exact overflow-threshold arithmetic quoted in the report |
| `raw/doc_excerpts.txt` | the two documented sentences, quoted from the pinned sources with locations |
| `logs/environment.txt` | tool versions, flags, archive hashes |
| `logs/commands.txt` | command lines and observed exit codes |
| `logs/*_build.log`, `logs/reference_*.log` | build and reference logs |
| `logs/*_raw_rerun.csv` | second witness runs, byte-identical to `raw/*raw.csv` |

## Pinned sources

Download into an empty `downloads/` directory and verify the hash before
extracting. These archives are untrusted data; do not run an interpreter or
build tool from inside the extracted trees.

| Source | URL | SHA-256 | Size |
| --- | --- | --- | --- |
| Boost 1.92.0 | `https://archives.boost.io/release/1.92.0/source/boost_1_92_0.tar.gz` | `c4a3b310ddd2472416e091067166b0713be97c63f38c212c484ada022fd296ce` | 235,384,471 |
| crlibm 1.0beta5, commit `eb3063791aa75bc9705b49283bf14250465220a7` | `https://github.com/taschini/crlibm/archive/eb3063791aa75bc9705b49283bf14250465220a7.tar.gz` | `79d626f9515fccdbc724e8b74cbf848331317b1891ef22aaaf6a195e20bbbfe5` | 2,432,556 |

Version stamps: `BOOST_VERSION 109200` and `BOOST_LIB_VERSION "1_92"`;
crlibm's `VERSION` file reads `1.0beta5`. The Boost hash also matches the
release site's own `.json` descriptor for 1.92.0.

## Requirements

- Windows with Visual Studio 18 Community (or another MSVC with the x64 C++
  tools). The build script calls
  `C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat`;
  change that one line if your installation is elsewhere.
- Python 3.10 or newer with mpmath 1.4.1 in the main `site-packages`. The
  scripts are run with `python -I`, which deliberately hides user
  `site-packages`; install mpmath into the interpreter's own environment if it
  is not already there. On Linux the interpreter is usually `python3`.
- WSL Ubuntu or a native Linux host with gcc, make, autoconf, automake and
  libtool. The rerun used gcc 15.2.0 on Ubuntu 26.04.

Boost is header-only, so there is no Boost build step.

## Reproduce

### 1. Fetch and verify the sources

From the root of this directory, in Git Bash or a POSIX shell:

```sh
mkdir -p downloads third_party raw logs
curl -L -o downloads/boost_1_92_0.tar.gz \
  https://archives.boost.io/release/1.92.0/source/boost_1_92_0.tar.gz
echo "c4a3b310ddd2472416e091067166b0713be97c63f38c212c484ada022fd296ce  downloads/boost_1_92_0.tar.gz" | sha256sum -c -
tar -xzf downloads/boost_1_92_0.tar.gz -C third_party

curl -L -o downloads/crlibm-eb306379.tar.gz \
  https://github.com/taschini/crlibm/archive/eb3063791aa75bc9705b49283bf14250465220a7.tar.gz
echo "79d626f9515fccdbc724e8b74cbf848331317b1891ef22aaaf6a195e20bbbfe5  downloads/crlibm-eb306379.tar.gz" | sha256sum -c -
tar -xzf downloads/crlibm-eb306379.tar.gz -C third_party
mv third_party/crlibm-eb3063791aa75bc9705b49283bf14250465220a7 \
   third_party/crlibm-1.0beta5
```

(`sha256sum` on macOS is `shasum -a 256`; the three hashes are the same.)

### 2. Boost.Math Airy

```sh
python -I scripts/make_boost_inputs.py
cmd //c scripts\\build_boost.bat
raw/boost_airy_witness.exe raw/boost_airy_inputs.txt > raw/boost_airy_raw.csv
python -I scripts/reference_boost.py
```

The input generator writes 2349 distinct arguments: the named points and
their immediate neighbours, dense windows of step 0.001 around each named
point, coarse windows around them, and the two audit grids reconstructed from
the row table in the deposit (`linear:65` over `[0.001, 100]` and `linear:33`
over `[30, 35.9]`). The reference script needs mpmath and writes
`raw/boost_airy_reference.csv` and `raw/boost_independent_check.txt`.

`cmd //c` is the Git Bash form; from `cmd.exe` run `scripts\build_boost.bat`
directly.

### 3. crlibm sinh

The library build and the witness run happen in WSL or Linux. The reference
script can run on either platform.

```sh
python -I scripts/make_crlibm_inputs.py        # may also be run in WSL
bash scripts/build_crlibm.sh                   # in WSL, from this directory
raw/crlibm_sinh_witness raw/crlibm_sinh_inputs.txt > raw/crlibm_sinh_raw.csv
python -I scripts/reference_crlibm.py
python -I scripts/crlibm_threshold_context.py
```

The build takes a few minutes. The driver reads 422 distinct arguments: the
named double, 64 neighbours on each side, the negative mirror, a scan of the
overflow tail, `DBL_MAX`, and a sanity point at 1.0. The reference script
writes `raw/crlibm_sinh_reference.csv` and
`raw/crlibm_sinh_independent_check.txt`; the threshold script writes
`raw/crlibm_threshold_context.txt`.

The configure flag `--disable-double-extended` matters. The audit's README
describes the non-double-extended configuration but its listed configure line
omits the flag, and `configure.ac` defaults the double-extended branch on. The
flag selects `log.c` rather than `log-de.c`; `csh_fast.c`, which contains
`sinh`, has no double-extended conditional.

### 4. Determinism

Repeat each witness run into `logs/` and compare byte for byte:

```sh
raw/boost_airy_witness.exe raw/boost_airy_inputs.txt > logs/boost_airy_raw_rerun.csv
cmp raw/boost_airy_raw.csv logs/boost_airy_raw_rerun.csv

raw/crlibm_sinh_witness raw/crlibm_sinh_inputs.txt > logs/crlibm_sinh_raw_rerun.csv
cmp raw/crlibm_sinh_raw.csv logs/crlibm_sinh_raw_rerun.csv
```

Both `cmp` commands exited 0 with no output at the rerun; the deposited
`logs/*_raw_rerun.csv` files are those second runs.

## Expected headline results

Boost.Math, error in units of eps = 2^-52, against mpmath at 60 digits:

| Function | Argument (double) | Bits | Error (eps) | Audit/paper |
| --- | --- | --- | --- | --- |
| `airy_ai` | 35.9 | `4041f33333333333` | 121.395340 | 121.395 |
| `airy_ai` | 78.1252 (literal) | `4053880346dc5d64` | 148.297779 | — |
| `airy_ai` | 78.12521875 (audit grid i=50) | `4053880395810625` | 414.878957 | — |
| `airy_bi` | 35.9 | `4041f33333333333` | 120.848885 | — |
| `airy_bi` | 78.1252 (literal) | `4053880346dc5d64` | 149.324772 | — |
| `airy_bi` | 78.12521875 (audit grid i=50) | `4053880395810625` | 414.110896 | 414.111 |
| `airy_bi` | 88.258434442270058 (gate context) | `4056108a309d6d6f` | 529.766809 | 529.767 |

The four-decimal label `78.1252` is the paper's rendering of the audit's
`linear:65` grid point `i = 50`, whose exact double is `78.12521875`. The
paper's `414.111` reproduces at that grid point. The literal double named by
the four-decimal label also breaches the sentence, at 149.325 eps.

crlibm at `0x408633ce8fb9f87e`:

| Mode | Result bits | Correct? | Exact value / required result |
| --- | --- | --- | --- |
| `sinh_rn` | `7ff0000000000000` (+inf) | yes | exact `1.7976931348623787561e308` is above the round-to-nearest overflow midpoint |
| `sinh_ru` | `7ff0000000000000` (+inf) | yes | positive overflow rounds up to +inf |
| `sinh_rd` | `7ff0000000000000` (+inf) | **no** | clause 7.4 requires `7fefffffffffffff` (+DBL_MAX) |
| `sinh_rz` | `7ff0000000000000` (+inf) | **no** | clause 7.4 requires `7fefffffffffffff` (+DBL_MAX) |

At the negative mirror `0xc08633ce8fb9f87e`, `sinh_rn` and `sinh_rd` are
correct and `sinh_ru` and `sinh_rz` wrongly return `fff0000000000000` (−inf)
where `ffefffffffffffff` (−DBL_MAX) is required. The remaining 420 sampled
inputs match in all four modes.

## What this deposit does not contain

`downloads/` (237 MB of upstream archives), `third_party/` (the extracted
Boost and crlibm trees, about 199 MB), and the built witness binaries and
object file are not included. They are third-party sources, they are large,
and the paper states that crlibm's source is not redistributed. The pins above
are the provenance: the hash fixes the bytes and the version stamps fix the
parsed source. Rebuilding is the first two commands of each reproduction
section.

`MANIFEST.sha256` lists every included file with its hash.

## Note on untrusted data

The archives and their extracted trees are data, not code to import. The
scripts live in `scripts/`, take paths as arguments, and are run with
`python -I`, which removes the script's directory and the current directory
from `sys.path`. Do not run Python from inside `third_party/` or `raw/`: a
planted `json.py` or similar would be imported by the interpreter.
