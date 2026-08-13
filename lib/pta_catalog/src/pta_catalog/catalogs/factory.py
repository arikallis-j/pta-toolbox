"""Factory functions for creating Catalog objects."""

from typing import Optional, Dict, List, Any 
from pathlib import Path

import numpy as np
import pandas as pd

from .catalog import Catalog
from .galactics import GALACTICS
from ..io.pickle_io import load_data, dump_data

def load_catalog(path: Path, name: str, prefix: bool = True) -> Catalog:
    """Load a catalog from storage."""
    stem = f"{name}_cat" if prefix else name
    data = load_data(path, stem)
    return Catalog(data, name=name)

def save_catalog(catalog: Catalog, path: Path, prefix: bool = True) -> None:
    """Save a catalog to storage."""
    stem = f"{catalog.name}_cat" if prefix else catalog.name
    dump_data(catalog.data, path, stem)

def make_catalog(n_psr: int, config: Dict[str, Any], mode: str = 'simple', name: Optional[str] = None) -> Catalog:
    name = name or f"{mode}-n{n_psr}"
    data = pd.DataFrame(np.nan, index=range(n_psr), columns=range(0))
    df = GALACTICS[mode](n_psr, config)
    data = pd.concat([data, df], axis=1)
    return Catalog(data, name=name)

def make_mixed_catalog(synth_catalog: Catalog, real_catalog: Catalog, seed: int = 42, fields: List[str] = ['F0','F1','PX','PMRA','PMDEC','PX','DM']):
    if len(synth_catalog) <= len(real_catalog):
        synth_df = synth_catalog.data
        real_df = real_catalog.sample(n_psr=len(synth_catalog), seed=seed).data
        name = f"{synth_catalog.name}_mixed_{real_catalog.name}-s{seed}"
    else:
        raise ValueError(f"len(synth_catalog) > len(real_catalog)")
    for col in fields:
        synth_df[col] = real_df[col].values
    synth_df['TYPE'] = 'MIX'
    return Catalog(data=synth_df, name=name)