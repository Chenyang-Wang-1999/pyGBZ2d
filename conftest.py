import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pytest


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: slow test (skipped by default, use --run-slow)")
