import glob
import numpy as np
import xarray as xr
from pathlib import Path
from enterprise.pulsar import Pulsar
from enterprise.signals.utils import createfourierdesignmatrix_red

def load_pulsars(path: Path) -> list:
    """Load all Pulsar objects found under path / {par,tim}."""
    pulsars = []
    parfiles = sorted(glob.glob(str(path / "par" / "*.par")))
    timfiles = sorted(glob.glob(str(path / "tim" / "*.tim")))
    for p, t in zip(parfiles, timfiles):
        pulsars.append(Pulsar(p, t, timing_package="tempo2"))
    return pulsars

def global_frequency_grid(pulsars: list, n_freqs: int = 30):
    """Compute the shared Tspan (max baseline) and frequency grid f_k = k / Tspan."""
    tspan = max(psr.toas.max() - psr.toas.min() for psr in pulsars)
    
    # Nyquist limit: conservative choice (worst pulsar sets the ceiling)
    dt_min = max(np.diff(np.sort(psr.toas)).min() for psr in pulsars)
    f_nyquist = 1.0 / (2.0 * dt_min)
    
    # Maximum number of modes within Nyquist
    n_freqs_max = int(f_nyquist * tspan) - 1
    if n_freqs > n_freqs_max:
        print(f"Warning: n_freqs={n_freqs} exceeds Nyquist limit {n_freqs_max}, clipping.")
        n_freqs = n_freqs_max
    elif n_freqs == -1:
        n_freqs = n_freqs_max
    
    freqs = np.arange(1, n_freqs + 1) / tspan
    return tspan, freqs, n_freqs

def project_to_fourier_domain(pulsars: list, n_freqs: int = 30):
    """
    Project all pulsars' residuals onto the shared Fourier basis.

    Returns
    -------
    freqs : ndarray, shape (n_freqs,)
        Frequency grid.
    r_freq : ndarray, shape (n_freqs, n_psr), complex128
        Complex Fourier amplitudes: r_tilde(f_k) = b_k - i*a_k
    N_cov  : ndarray, shape (n_freqs, n_psr, n_psr), complex128
        Noise covariance matrix per frequency bin.
        Diagonal: per-pulsar noise variance from GLS.
        Off-diagonal: zero (no inter-pulsar correlations in noise).
    """
    
    tspan, freqs, n_freqs = global_frequency_grid(pulsars, n_freqs=n_freqs)

    n_psr = len(pulsars)
    r_freq = np.zeros((n_freqs, n_psr), dtype=complex)
    N_cov  = np.zeros((n_freqs, n_psr, n_psr), dtype=complex)

    T = np.array([[-1j, 1.0]])  # shape (1, 2), converts [a_k, b_k] -> b_k - i*a_k

    for i, psr in enumerate(pulsars):
        F, _ = createfourierdesignmatrix_red(psr.toas, nmodes=n_freqs, Tspan=tspan)

        weights  = 1.0 / psr.toaerrs**2
        FtCinv   = F.T * weights
        FtCinvF  = FtCinv @ F
        FtCinvr  = FtCinv @ psr.residuals

        N_real = np.linalg.inv(FtCinvF)  # (2*n_freqs, 2*n_freqs)
        r_real = N_real @ FtCinvr        # (2*n_freqs,)

        for k in range(n_freqs):
            ab    = r_real[2*k : 2*k+2]
            block = N_real[2*k:2*k+2, 2*k:2*k+2]

            r_freq[k, i]       = (T @ ab)[0]
            N_cov[k, i, i]     = (T @ block @ T.conj().T)[0, 0]  # diagonal block only

    return freqs, r_freq, N_cov

def load_pulsar_coords(pulsars: list):
    theta = np.array([psr.theta for psr in pulsars])
    phi = np.array([psr.phi   for psr in pulsars])
    return theta, phi

def load_pulsar_names(pulsars: list):
    names = np.array([psr.name for psr in pulsars])
    return names

def make_observations(pulsars: list, name: str = "sample", n_freqs: int = -1):
    tspan, freqs, n_freqs = global_frequency_grid(pulsars, n_freqs=n_freqs)
    freqs, r_freq, N_cov = project_to_fourier_domain(pulsars, n_freqs=n_freqs)
    theta, phi = load_pulsar_coords(pulsars)
    names = load_pulsar_names(pulsars)

    ds = xr.Dataset(
        {
            "r_freq": (["freq", "psr"], r_freq),
            "N_cov":  (["freq", "psr_i", "psr_j"], N_cov),
        },
        coords={
            "freq":   freqs,
            "psr":    names,
            "psr_i":    names,
            "psr_j":    names,
            "theta": (["psr"], theta),  
            "phi":   (["psr"], phi),   
        },
        attrs={
            "name":  name,
            "tspan":  tspan,
            "n_freqs": n_freqs,
        }
    )
    
    return ds

