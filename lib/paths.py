"""Paths shared by every Python script, from the same environment variables as config.sh."""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RUN = Path(os.environ.get("GI_RUN", Path.home() / "gbm_imprinting_run"))
DATA = Path(os.environ.get("GI_DATA", RUN / "data"))           # downloaded inputs
WORK = Path(os.environ.get("GI_WORK", RUN / "work"))           # intermediate files
RESULTS = Path(os.environ.get("GI_RESULTS", RUN / "results"))  # final tables and figures
SHIPPED = REPO / "data"                                        # small inputs in the repository
FIG = RESULTS / "figures"
TAB = RESULTS / "tables"
