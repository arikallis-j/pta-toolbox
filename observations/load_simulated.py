"""Load pulsars and project their residuals onto a shared Fourier basis.

Mirrors the pulsar-loading pattern in observations/main.py, but instead of
working in the time domain, this produces per-pulsar frequency-domain
residuals and noise covariances on a global frequency grid, ready to feed
into the SVD sky-mapping pipeline in scripts/main.py.
"""
import glob

import numpy as np
import typer

import ptatoolbox as pta
from enterprise.pulsar import Pulsar
from enterprise.signals.utils import createfourierdesignmatrix_red

SAMPLES = ("psr_test_1", "psr_test_2")


def load_pulsars(dm: pta.utils.DataManager, sample: str) -> list:
    """Load all Pulsar objects found under dm.raw / <sample> / {par,tim}."""
    pulsars = []
    datapath = dm.raw / sample
    parfiles = sorted(glob.glob(str(datapath / "par" / "*.par")))
    timfiles = sorted(glob.glob(str(datapath / "tim" / "*.tim")))
    for p, t in zip(parfiles, timfiles):
        pulsars.append(Pulsar(p, t, timing_package="tempo2"))
    return pulsars


def global_frequency_grid(pulsars: list, n_freqs: int = 30):
    """Compute the shared Tspan (max baseline) and frequency grid f_k = k / Tspan."""
    tspan = max(psr.toas.max() - psr.toas.min() for psr in pulsars)
    freqs = np.arange(1, n_freqs + 1) / tspan
    return tspan, freqs


def project_to_fourier_domain(psr: Pulsar, tspan: float, n_freqs: int):
    """GLS-project one pulsar's time-domain residuals onto the shared Fourier basis.

    Fourier design matrix columns alternate sin/cos: column 2k is sin(f_k) and
    column 2k + 1 is cos(f_k), for k = 0, ..., n_freqs - 1 (enterprise's
    convention). The complex spectrum can be recovered from the returned
    r_freq as r_tilde(f_k) = r_freq[2k + 1] - 1j * r_freq[2k].

    Returns
    -------
    r_freq : ndarray, shape (2 * n_freqs,)
        (F^T C^-1 F)^-1 F^T C^-1 r_time
    N_cov : ndarray, shape (2 * n_freqs, 2 * n_freqs)
        (F^T C^-1 F)^-1
    """
    F, _ = createfourierdesignmatrix_red(psr.toas, nmodes=n_freqs, Tspan=tspan)

    # TODO: project out the timing model here (G = I - M(M^T M)^-1 M^T applied
    # to F and r_time via psr.Mmat) before the GLS fit below.

    weights = 1.0 / psr.toaerrs**2
    FtCinv = F.T * weights
    FtCinvF = FtCinv @ F

    N_cov = np.linalg.inv(FtCinvF)
    r_freq = N_cov @ (FtCinv @ psr.residuals)
    return r_freq, N_cov


def load_simulated_data(sample: str, n_freqs: int = 30, root_dir: str = "./data") -> dict:
    """Load pulsars from a single `sample` and build frequency-domain data on a shared grid."""
    if sample not in SAMPLES:
        raise ValueError(f"unknown sample {sample!r}, expected one of {SAMPLES}")

    dm = pta.utils.DataManager(root_dir)
    pulsars = load_pulsars(dm, sample)

    tspan, freqs = global_frequency_grid(pulsars, n_freqs)

    positions = np.array([psr.pos for psr in pulsars])
    r_freq = np.zeros((len(pulsars), 2 * n_freqs))
    N_cov = np.zeros((len(pulsars), 2 * n_freqs, 2 * n_freqs))

    for i, psr in enumerate(pulsars):
        r_freq[i], N_cov[i] = project_to_fourier_domain(psr, tspan, n_freqs)

    return {
        "names": [psr.name for psr in pulsars],
        "freqs": freqs,
        "positions": positions,
        "r_freq": r_freq,
        "N_cov": N_cov,
    }


def main(sample: str = "psr_test_1", n_freqs: int = 30):
    data = load_simulated_data(sample, n_freqs=n_freqs)

    n_psr = data["positions"].shape[0]
    freqs = data["freqs"]
    print(f"sample: {sample}")
    print(f"pulsars loaded: {n_psr}")
    print(f"frequency grid: {n_freqs} modes, {freqs[0]:.3e} - {freqs[-1]:.3e} Hz")
    print(f"positions: {data['positions'].shape}")
    print(f"r_freq:    {data['r_freq'].shape}")
    print(f"N_cov:     {data['N_cov'].shape}")

if __name__ == "__main__":
    typer.run(main)
