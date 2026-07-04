"""Factory functions for creating Catalog objects."""

from typing import Optional, Dict, List, Any 
from pathlib import Path

from .catalog import Catalog
from .models import make_synthetics
from ..io.manager import DataManager
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

def make_synthetic_catalog(n_psr: int, name: Optional[str] = None, method: Optional[str] = 'test', params: Optional[Dict[str, Any]] = None) -> Catalog:
    """Create a fully synthetic catalog."""
    params = params or {}
    name = name or 'sample'
    data = make_synthetics(n_psr, method, params)
    return Catalog(data, name=name)

def make_mixed_catalog(
    real_catalog: Catalog, 
    n_psr: Optional[int] = None,
    name: Optional[str] = None,
    method: str = 'test', 
    params: Optional[Dict[str, Any]] = None, 
    seed: int = 42, 
    fields: Optional[List[str]] = ['PSRJ', 'RAJD', 'DECJD'], 
) -> Catalog:
    """Create a mixed catalog by replacing specified fields with synthetic values."""
    if name is None:
        name =  f"mixed_{real_catalog.name}"
    if n_psr is None:
        n_psr = len(real_catalog)
        sampled_df = real_catalog.data
        final_name = f"{name}"
    else:
        sampled_df = real_catalog.sample(n_psr=n_psr, seed=seed).data
        final_name = f"{name}-n{n_psr}-s{seed}"
    
    synth_df = make_synthetics(n_psr, method, params)

    for col in fields:
        if col in synth_df.columns:
            sampled_df[col] = synth_df[col].values
        else:
            sampled_df[col] = float('nan')
    
    return Catalog(data=sampled_df, name=final_name)


def make_catalog(
    n_psr: int,
    real_catalog: Optional[Catalog] = None, 
    name: Optional[str] = None,
    method: str = 'test', 
    params: Optional[Dict[str, Any]] = None, 
    seed: int = 42, 
    fields: Optional[List[str]] = ['PSRJ', 'RAJD', 'DECJD'], 
) -> Catalog:
    """Create a catalog."""
    if real_catalog is None:
        catalog = make_synthetic_catalog(
            n_psr=n_psr, 
            name=name, 
            method=method, 
            params=params
        )
    else:
        catalog = make_mixed_catalog(
            real_catalog=real_catalog,
            n_psr=n_psr,
            name=name,
            method=method,
            params=params,
            seed=seed,
            fields=fields,
        )

    return catalog