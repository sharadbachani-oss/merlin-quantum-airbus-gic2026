# -*- coding: utf-8 -*-
"""Airbus Dirac closed-orbit even-parity leftover. Jennings NS bound stays."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import run_dirac_advantage as R
if __name__ == "__main__":
    raise SystemExit(R.main())
