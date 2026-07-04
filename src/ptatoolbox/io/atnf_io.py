"""I/O operations for ATNF catalog."""

import psrqpy
import pandas as pd
from typing import Optional
from pathlib import Path

from .pickle_io import load_data, dump_data
from .manager import DataManager
from ..core.constants import ATNF_FORMAT

ATNF_STEM = "atnf"
CUT_ATNF_STEM = "cut_atnf"

def download_atnf(path: Path) -> pd.DataFrame:
    """Download ATNF catalog and cache it in storage."""
    atnf = psrqpy.QueryATNF().pandas
    dump_data(atnf, path, ATNF_STEM)
    return atnf
    
def normalize_atnf_data(data: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Ensure DataFrame has all required columns from ATNF_FORMAT."""
    required_cols = list(ATNF_FORMAT.values())
    if data is None:
        return pd.DataFrame(columns=required_cols)   
    df = data.copy()
    for col in required_cols:
        if col not in df.columns:
            if col == 'PSRJ':
                df[col] = [f"S{i}" for i in range(len(df))]
            else:
                df[col] = float('nan')
    return df[required_cols]