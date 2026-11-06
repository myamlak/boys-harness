"""Resolve the paths the certified-bound record reads, and the run directory
it writes.

Every path is taken from the command line when given, else from the
documented environment variable, else from a default relative to this file:

    kernel checkout     --kernel        BOYS_ROOT         <repo>/../boys
    sweep data          --sweep         CERT_SWEEP        ../amplification-form
    model record        --model-record  CERT_MODEL_RECORD inputs/model-record
    exponential bits    --exp-bits      CERT_EXP_BITS     exp/exp_bits.txt
    deposit             --deposit       CERT_DEPOSIT      ../apriori-bound
    run directory       --out           -                 (required)

The defaults assume the layout of the harness checkout; the kernel default
is the CI convention of a boys clone beside the harness.  Nothing is written
outside the directory selected with --out.

The deposit (amplification-sweeps/apriori-bound) imports its own siblings
(``import bounds``, ``import coeffs``), so ``deposit_on_path`` puts that one
in-repo directory on sys.path for the duration of the import.  No external
working copy is ever added to sys.path.
"""
import contextlib
import os
import sys
from pathlib import Path

RECORD_DIR = Path(__file__).resolve().parent.parent
CODE_DIR = RECORD_DIR / "code"
EXPLORE_DIR = RECORD_DIR / "explore"
REPO_ROOT = RECORD_DIR.parent.parent

DEFAULT_KERNEL = REPO_ROOT.parent / "boys"
DEFAULT_SWEEP = RECORD_DIR.parent / "amplification-form"
DEFAULT_MODEL_RECORD = RECORD_DIR / "inputs" / "model-record"
DEFAULT_EXP_BITS = RECORD_DIR / "exp" / "exp_bits.txt"
DEFAULT_DEPOSIT = RECORD_DIR.parent / "apriori-bound"

_KERNEL_HINT = (
    "pass --kernel PATH or set BOYS_ROOT to a boys checkout at v5.2.5 "
    "(git clone --branch v5.2.5 https://github.com/myamlak/boys.git ../boys)")


def _select(value, env, default):
    if value is not None:
        return Path(value).expanduser().resolve()
    env_value = os.environ.get(env)
    if env_value:
        return Path(env_value).expanduser().resolve()
    return Path(default).resolve()


def kernel(value=None):
    path = _select(value, "BOYS_ROOT", DEFAULT_KERNEL)
    if not (path / "include" / "boys" / "boys_coefficients.hpp").is_file():
        raise SystemExit(f"kernel checkout not found at {path}\n{_KERNEL_HINT}")
    return path


def sweep(value=None):
    path = _select(value, "CERT_SWEEP", DEFAULT_SWEEP)
    needed = ("data/reference.csv", "results/cells.csv", "results/per_order.csv")
    missing = [name for name in needed if not (path / name).is_file()]
    if missing:
        raise SystemExit(
            f"sweep data not found at {path} (missing {', '.join(missing)})\n"
            "pass --sweep PATH or set CERT_SWEEP to a sweep root with data/ "
            "and results/")
    return path


def model_record(value=None):
    path = _select(value, "CERT_MODEL_RECORD", DEFAULT_MODEL_RECORD)
    for name in ("prediction.csv", "prediction_record.csv"):
        if not (path / name).is_file():
            raise SystemExit(
                f"model record not found at {path}/{name}\n"
                "pass --model-record DIR or set CERT_MODEL_RECORD to a "
                "directory holding prediction.csv and prediction_record.csv")
    return path


def exp_bits(value=None):
    path = _select(value, "CERT_EXP_BITS", DEFAULT_EXP_BITS)
    if not path.is_file():
        raise SystemExit(
            f"exponential bit file not found at {path}\n"
            "pass --exp-bits FILE or set CERT_EXP_BITS")
    return path


def deposit(value=None):
    path = _select(value, "CERT_DEPOSIT", DEFAULT_DEPOSIT)
    if not (path / "bounds.py").is_file():
        raise SystemExit(
            f"a priori-bound deposit not found at {path}\n"
            "pass --deposit DIR or set CERT_DEPOSIT")
    return path


def out_dir(value):
    if value is None:
        raise SystemExit(
            "--out is required: the run writes only beneath the selected "
            "directory")
    path = Path(value).expanduser().resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path


def add_kernel(parser):
    parser.add_argument(
        "--kernel", metavar="DIR",
        help="boys checkout at v5.2.5 (default: $BOYS_ROOT, else <repo>/../boys)")


def add_inputs(parser, model_record_arg=False, deposit_arg=False,
               exp_bits_arg=True):
    parser.add_argument(
        "--sweep", metavar="DIR",
        help="sweep root with data/reference.csv and results/ "
             "(default: $CERT_SWEEP, else ../amplification-form)")
    if model_record_arg:
        parser.add_argument(
            "--model-record", metavar="DIR",
            help="directory with prediction.csv and prediction_record.csv "
                 "(default: $CERT_MODEL_RECORD, else inputs/model-record)")
    if deposit_arg:
        parser.add_argument(
            "--deposit", metavar="DIR",
            help="a priori-bound deposit directory "
                 "(default: $CERT_DEPOSIT, else ../apriori-bound)")
    if exp_bits_arg:
        parser.add_argument(
            "--exp-bits", metavar="FILE",
            help="exponential probe output (default: $CERT_EXP_BITS, else "
                 "exp/exp_bits.txt)")


def add_out(parser):
    parser.add_argument(
        "--out", metavar="DIR", required=True,
        help="run directory; every output is written beneath it")


@contextlib.contextmanager
def deposit_on_path(path):
    """Import the in-repo deposit's sibling-importing modules from its own
    directory.  ``path`` is resolved from the record or given on the command
    line; it is never an external working copy."""
    entry = str(Path(path).resolve())
    sys.path.insert(0, entry)
    try:
        yield
    finally:
        try:
            sys.path.remove(entry)
        except ValueError:
            pass
