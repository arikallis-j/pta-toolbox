"""Catalog class and related factories."""

import pandas as pd
from typing import Optional, List, Dict, Any

from .pulsar import Pulsar, make_pulsars, make_data
from ..core.constants import TEMPO_FORMAT
from ..io.atnf_io import normalize_atnf_data

class Catalog:
    """A container for pulsar catalog data with standardized columns."""
    def __init__(self, data: Optional[pd.DataFrame] = None, name: str = 'test'):
        """Initialize a Catalog."""
        self.name = name
        self.data = normalize_atnf_data(data)
    
    def pulsars(self) -> List[Pulsar]:
        """Convert catalog data to a list of Pulsar objects."""
        return make_pulsars(self.data)
    
    def tempo(self) -> pd.DataFrame:
        """Return data with columns renamed for TEMPO format."""
        return self.data.rename(columns=TEMPO_FORMAT)
    
    def add(self, pulsars: List[Pulsar]) -> 'Catalog':
        """Add pulsars to the catalog and return a new Catalog instance."""
        new_data = pd.concat([self.data, make_data(pulsars)], ignore_index=True)
        return Catalog(data=new_data, name=self.name)
    
    def sample(self, n_psr: int, seed: int = 42, replace: bool = False) -> 'Catalog':
        """Return a new catalog with a random sample of n pulsars."""
        sample_data = self.data.sample(n_psr, random_state=seed, replace=replace)
        new_name = f"{self.name}-n{n_psr}-s{seed}" + ("-rep" if replace else "")
        return Catalog(data=sample_data, name=new_name)
    
    def filter_by_names(self, names: List[str]) -> 'Catalog':
        """Return a new catalog containing only pulsars whose names are in the list."""
        mask = self.data['PSRJ'].isin(names)
        return Catalog(data=self.data[mask].copy(), name=f"{self.name}_filtered")
    
    def copy(self) -> 'Catalog':
        """Return a deep copy of the catalog."""
        return Catalog(data=self.data.copy(), name=self.name)
    
    def __repr__(self) -> str:
        return f"{self.name} catalog\n" + repr(self.data)
    
    def __str__(self) -> str:
        return f"{self.name} catalog\n" + str(self.data)

    def __len__(self) -> int:
        return len(self.data)