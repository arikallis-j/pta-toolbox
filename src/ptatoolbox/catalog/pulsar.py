"""Pulsar class and related factories."""

import pandas as pd

from typing import NamedTuple, List
from ..core.constants import ATNF_FORMAT
from ..io.atnf_io import normalize_atnf_data

class Pulsar(NamedTuple):
    """Immutable container for pulsar parameters."""
    name: str
    ra: float                   # deg
    dec: float                  # deg
    f0: float = 1.0             # Hz
    f1: float = float('nan')    # s^-2
    pmra: float = float('nan')  # mas yr^-1
    pmdec: float = float('nan') # mas yr^-1
    px: float = float('nan')    # mas
    dm: float = float('nan')    # cm^-3 pc

def make_pulsar(params: dict) -> Pulsar:
    """Create a Pulsar object from a dictionary of parameters."""
    return Pulsar(**params)

def make_pulsars(data: pd.DataFrame) -> List[Pulsar]:
    """Convert a DataFrame (with ATNF column names) to a list of Pulsar objects."""
    subset = normalize_atnf_data(data)
    subset.rename(columns={v: k for k, v in ATNF_FORMAT.items()}, inplace=True)
    return [Pulsar(**row) for row in subset.to_dict('records')]

def make_data(pulsars: List[Pulsar]) -> pd.DataFrame:
    """Convert a list of Pulsar objects to a DataFrame with ATNF columns."""
    df = pd.DataFrame(pulsars)
    df.rename(columns=ATNF_FORMAT, inplace=True)
    return df
