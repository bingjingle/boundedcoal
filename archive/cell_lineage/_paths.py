"""Paths for the archived cell-lineage analysis."""
import os
from pathlib import Path

ARCHIVE_ROOT = Path(__file__).resolve().parent
DATA_DIR = Path(os.environ.get("BOUNDEDCOAL_CELL_LINEAGE_DATA_DIR", ARCHIVE_ROOT)).expanduser()
OUTPUT_DIR = Path(os.environ.get("BOUNDEDCOAL_CELL_LINEAGE_OUTPUT_DIR", ARCHIVE_ROOT.parents[1] / "outputs" / "cell_lineage")).expanduser()


def data_file(name):
    return str(DATA_DIR / name)


def output_file(name):
    path = OUTPUT_DIR / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return str(path)
