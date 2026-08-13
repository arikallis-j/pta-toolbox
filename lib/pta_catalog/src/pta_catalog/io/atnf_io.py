"""I/O operations for ATNF catalog."""

import psrqpy
import numpy as np
import pandas as pd
from typing import Optional
from pathlib import Path
from collections import Counter


from .pickle_io import load_data, dump_data
from ..core.constants import ATNF_FORMAT

ATNF_STEM = "atnf"
CUT_ATNF_STEM = "cut_atnf"

def download_atnf(path: Path) -> pd.DataFrame:
    """Download ATNF catalog and cache it in storage."""
    atnf = psrqpy.QueryATNF().pandas
    dump_data(atnf, path, ATNF_STEM)
    return atnf
    
def normalize_atnf_data(data: Optional[pd.DataFrame] = None, prefix: Optional[str] = None) -> pd.DataFrame:
    """Ensure DataFrame has all required columns from ATNF_FORMAT."""
    required_cols = list(ATNF_FORMAT.values())
    if data is None:
        return pd.DataFrame(columns=required_cols)   
    df = data.copy()
    for col in required_cols:
        if col not in data.columns:
            if col == 'PSRJ':
                if ('RAJD' not in data.columns) and ('DECJD' not in data.columns):
                    df[col] = get_simple_names(range(len(df)), prefix=prefix)
                else:
                    df[col] = get_atnf_names(data['RAJD'], data['DECJD'], prefix=prefix)
            elif col == 'TYPE':
                df[col] = 'REL'
            else:
                df[col] = float('nan')
    return df[required_cols]

def get_simple_name(idx: int, prefix: Optional[str] = None):
    prefix = 'J' or prefix
    name = f"{prefix}{idx}"
    return name

def get_atnf_name(ra: float, dec: float, prefix: Optional[str] = None):
    prefix = 'J' or prefix
    ra_h = np.floor(ra/15.0).astype(int)
    ra_m = np.floor(np.round((ra/15.0 - ra_h) * 60, decimals=1)).astype(int)
    dec_d = np.floor(np.abs(dec)).astype(int)
    dec_m = np.floor(np.round((np.abs(dec) - dec_d) * 60, decimals=1)).astype(int)
    dec_s = "+" if np.sign(dec)>=0.0 else "-"
    name = f"{prefix}{ra_h:02d}{ra_m:02d}{dec_s}{dec_d:02d}{dec_m:02d}"
    return name

def get_simple_names(idx: int, prefix: Optional[str] = None):
    base_names = [get_simple_name(k, prefix) for k in idx]
    return base_names

def get_atnf_names(ra: pd.DataFrame, dec: pd.DataFrame, prefix: Optional[str] = None):
    base_names = [get_atnf_name(r, d, prefix) for r, d in zip(ra, dec)]
    freq = Counter(base_names)
    counters = {name: 0 for name in freq}
    final_names = []
    for name in base_names:
        cnt = counters[name]
        if cnt == 0:    # 'coords': {
    #     'model': 'sphere',
    #     'params': {}
    # },
            final_names.append(name)
        else:
            if cnt <= 26:
                suffix = chr(ord('A') + cnt - 1)
            else:
                idx = cnt - 27
                first = chr(ord('a') + idx // 26)
                second = chr(ord('a') + idx % 26)
                suffix = first + second
            final_names.append(name + suffix)
        counters[name] += 1
    
    return final_names