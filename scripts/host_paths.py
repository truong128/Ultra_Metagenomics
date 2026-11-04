"""Machine-specific locations.

The repository does not store reference databases. Set METAX_DB_ROOT when the
database directory is not the `databases` folder next to this repository.
Set METAX_CONDA when Miniforge is not in the home directory.
"""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def database_root():
    env = os.environ.get('METAX_DB_ROOT', '').strip()
    if env:
        return Path(env).expanduser()
    sibling = REPO.parent / 'databases'
    if sibling.is_dir():
        return sibling
    return REPO / 'databases'


def conda_root():
    env = os.environ.get('METAX_CONDA', '').strip()
    if env:
        return Path(env).expanduser()
    return Path.home() / 'miniforge3'
