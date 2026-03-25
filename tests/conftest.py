"""Pytest fixtures for FusionSQL tests."""

import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def sample_ariba_path():
    """Path to the sample Ariba fusion TSV fixture."""
    return Path(__file__).parent / "fixtures" / "ariba_fusions.tsv"


@pytest.fixture
def sample_starfusion_path():
    """Path to the sample StarFusion fusion TSV fixture."""
    return Path(__file__).parent / "fixtures" / "starfusion_fusions.tsv"


@pytest.fixture
def temp_output_dir(tmp_path):
    """Temporary directory for output files."""
    return tmp_path / "output"


@pytest.fixture
def temp_output_dir_created(temp_output_dir):
    """Creates the temp output directory."""
    temp_output_dir.mkdir(parents=True, exist_ok=True)
    return temp_output_dir
