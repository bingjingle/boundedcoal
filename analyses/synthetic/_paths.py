"""Local paths shared by the historical synthetic experiment scripts."""
import os
from pathlib import Path

SYNTHETIC_ROOT = Path(__file__).resolve().parent
DATA_DIR = Path(os.environ.get("BOUNDEDCOAL_SYNTHETIC_DATA_DIR", SYNTHETIC_ROOT / "data")).expanduser()
OUTPUT_DIR = Path(os.environ.get("BOUNDEDCOAL_SYNTHETIC_OUTPUT_DIR", SYNTHETIC_ROOT.parents[1] / "outputs" / "synthetic")).expanduser()


def data_file(name):
    return str(DATA_DIR / name)


def output_file(name):
    path = OUTPUT_DIR / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return str(path)
