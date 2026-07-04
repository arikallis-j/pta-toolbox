"""Low-level I/O operations for pickle files."""

import pandas as pd
from pathlib import Path
from typing import Union

PICKLE_EXT = ".pkl"

def load_data(path: Union[str, Path], stem: str) -> pd.DataFrame:
    """Load a pandas DataFrame from a pickle file."""
    path = Path(path)
    return pd.read_pickle(path / f"{stem}{PICKLE_EXT}")

def dump_data(data: pd.DataFrame, path: Union[str, Path], stem: str) -> None:
    """Save a pandas DataFrame to a pickle file."""
    path = Path(path)
    data.to_pickle(path / f"{stem}{PICKLE_EXT}")