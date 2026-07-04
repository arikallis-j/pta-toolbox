"""Synthetic pulsar catalog generators."""

import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional

from .funcs import (
    psr2sph, sph2psr, get_names,
    isotropic_sphere, isotropic_ball, 
    isotropic_cap, isotropic_cone, isotropic_ring,
)
from .pulsar import Pulsar, make_data

def _make_pulsar_dicts(Phi: np.ndarray, Theta: np.ndarray, Rho: np.ndarray, prefix: str = 'S') -> List[Dict[str, Any]]:
    """Convert spherical coordinates to list of pulsar parameter dicts."""
    Ra, Dec, Px = sph2psr(Phi, Theta, Rho)
    names = get_names(Ra, Dec, prefix=prefix)
    return [
        {'name': name, 'ra': ra, 'dec': dec, 'px': px}
        for name, ra, dec, px in zip(names, Ra, Dec, Px)
    ]

def simple_test(n_psr: int) -> List[Dict[str, Any]]:
    """Simplest generator: isotropic sphere with default parameters, no px."""
    Phi, Theta, Rho = isotropic_sphere(n_psr, 42, np.nan)
    Ra, Dec, _ = sph2psr(Phi, Theta, Rho)
    names = [f"S{k}" for k in range(n_psr)]
    return [
        {'name': name, 'ra': ra, 'dec': dec, 'px': np.nan}
        for name, ra, dec in zip(names, Ra, Dec)
    ]

def simple_sphere(n_psr: int, seed_psr: int = 42, radius: float = np.nan) -> List[Dict[str, Any]]:
    """Isotropic distribution on a sphere."""
    Phi, Theta, Rho = isotropic_sphere(n_psr, seed_psr, radius)
    return _make_pulsar_dicts(Phi, Theta, Rho)

def simple_ball(n_psr: int, seed_psr: int = 42, radius: float = 1.0) -> List[Dict[str, Any]]:
    """Uniform distribution inside a ball."""
    Phi, Theta, Rho = isotropic_ball(n_psr, seed_psr, radius)
    return _make_pulsar_dicts(Phi, Theta, Rho)

def simple_cone(n_psr: int, seed_psr: int = 42, alpha: float = 10.0, ra_0: float = 0.0, dec_0: float = 0.0, radius: float = 1.0) -> List[Dict[str, Any]]:
    """Distribution within a cone."""
    phi_0, theta_0, _ = psr2sph(ra_0, dec_0)
    Phi, Theta, Rho = isotropic_cone(n_psr, seed_psr, radius, phi_0, theta_0, np.deg2rad(alpha))
    return _make_pulsar_dicts(Phi, Theta, Rho)

def simple_cap(n_psr: int, seed_psr: int = 42, alpha: float = 10.0, ra_0: float = 0.0, dec_0: float = 0.0, radius: float = np.nan) -> List[Dict[str, Any]]:
    """Distribution on a spherical cap."""
    phi_0, theta_0, _ = psr2sph(ra_0, dec_0)
    Phi, Theta, Rho = isotropic_cap(n_psr, seed_psr, radius, phi_0, theta_0, np.deg2rad(alpha))
    return _make_pulsar_dicts(Phi, Theta, Rho)

def simple_ring(n_psr: int, seed_psr: int = 42, alpha: float = 10.0, ra_0: float = 0.0, dec_0: float = 0.0, radius: float = np.nan) -> List[Dict[str, Any]]:
    """Distribution on a ring (circle) on the sphere."""
    phi_0, theta_0, _ = psr2sph(ra_0, dec_0)
    Phi, Theta, Rho = isotropic_ring(n_psr, seed_psr, radius, phi_0, theta_0, np.deg2rad(alpha))
    return _make_pulsar_dicts(Phi, Theta, Rho)

def make_synthetics(n_psr: int, method: str, params: Optional[Dict[str, Any]] = None) -> pd.DataFrame:
    """Generate a synthetic catalog as a DataFrame with ATNF columns."""
    params = params or {}
    synth_data = methods[method](n_psr, **params)
    pulsars = [Pulsar(**item) for item in synth_data]
    return make_data(pulsars)

# Registry of available synthetic generators
methods = {
    'test': simple_test,
    'sphere': simple_sphere,
    'ball': simple_ball,
    'cap': simple_cap,
    'ring': simple_ring,
    'cone': simple_cone,
}