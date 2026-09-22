"""Frequency-resolved GW sky reconstruction from real pulsar data.

This replays the SVD dirty-map / Fisher-matrix / clean-map chain used in
`scripts/main.py`'s synthetic toy model, but on real per-pulsar frequency-domain
residuals from `observations/load_simulated.py`, at each Fourier frequency
independently.

Conventions (verified against the current repo, not assumed):

- **Spherical harmonics.** `scripts/main.py` uses complex harmonics via
  `scipy.special.sph_harm_y(l, m, theta, phi)` (theta = colatitude, phi =
  azimuth) directly, not a real-harmonic (+m/-m combined) basis. We replicate
  that exactly. `positions` (unit vectors from `load_simulated.py`) are
  converted to (theta, phi) with `healpy.vec2ang` (there is no
  `healpy.sphere2vec`).

- **Nl weighting / a_lm vs psi basis.** `ptatoolbox.utils.transform.calc_Nl`
  gives Nl(l) = sqrt(2 / ((l-1) l (l+1) (l+2))) = sqrt(2 (l-2)! / (l+2)!) for
  l > 1. `scripts/main.py` generates its synthetic truth directly as
  `r_true = D_raw @ psi_vec` (D_raw = raw, unweighted 4*pi*Y_lm matrix acting
  on the *potential* psi_lm), then fits using the reweighted matrix
  `D = D_raw * (Nl/2)`. Solving `D @ x = r` for x therefore recovers
  `x = (2/Nl) * psi = alm` (matches `plm2alm`/`alm2plm` in transform.py, and
  matches `scripts/main.py` naming `a_vec_est` -> `alm_est` -> `plm_est =
  alm_est * Nl/2`). So the SVD/Fisher solve's *direct* output is the
  a_lm-normalized coefficient vector, and psi_lm = alm * Nl / 2 is derived
  from it — the reverse of `alm = 2*psi_lm/Nl`. We follow that (verified)
  direction here rather than deriving alm from psi_lm, since the latter would
  invert the actual relationship encoded by `scripts/main.py`'s generation
  step and by `transform.alm2plm`/`plm2alm`.

- **Fourier sign convention.** `scripts/main.py` has no time-domain Fourier
  step at all (it works with a single synthetic complex snapshot per
  pulsar, no `e^{i2*pi*f*t}` anywhere). The sign convention used here instead
  follows `observations/load_simulated.py`'s documented convention: Fourier
  design-matrix columns alternate sin (even) / cos (odd) per
  `enterprise.signals.utils.createfourierdesignmatrix_red`, and the complex
  residual is reconstructed as `r_tilde(f_k) = r_cos_k - 1j * r_sin_k`
  (r(t) = Re[r_tilde * e^{i*2*pi*f*t}]).

- **A(f) factor.** Neither `scripts/main.py` nor anything else in this repo
  implements `A(f) = 1/(i*2*pi*f)` — the toy model is frequency-independent.
  It is included here as a separate scalar factor multiplying the raw
  (unweighted) design matrix, per the physical relation `r_tilde(f) = A(f) *
  D_raw(f) @ psi_lm` relating actual frequency-domain residuals to the sky
  potential, before the Nl/2 reweighting is applied for the SVD/Fisher solve
  exactly as `scripts/main.py` reweights its (frequency-independent) D_raw.
"""
from pathlib import Path
from typing import Optional

import healpy as hp
import numpy as np
import typer
from scipy.special import sph_harm_y as sph_harm

from observations.load_simulated import load_simulated_data

L_MAX_DEFAULT = 8


def build_lm_index(l_max: int):
    """List of (l, m) pairs for l = 2..l_max, m = -l..l (scripts/main.py's l>=2 cut)."""
    return [(l, m) for l in range(2, l_max + 1) for m in range(-l, l + 1)]


def calc_Nl_vec(lm_index):
    """Nl(l) = sqrt(2 / ((l-1) l (l+1) (l+2))), broadcast over the (l, m) index list."""
    l_arr = np.array([l for l, _ in lm_index], dtype=np.float64)
    return np.sqrt(2.0 / ((l_arr - 1) * l_arr * (l_arr + 1) * (l_arr + 2)))


def build_raw_design_matrix(theta, phi, lm_index):
    """D_raw[i, nu] = 4*pi * Y_lm(theta_i, phi_i), complex Y_lm via scipy sph_harm_y."""
    l_arr = np.array([l for l, _ in lm_index])[None, :]
    m_arr = np.array([m for _, m in lm_index])[None, :]
    Y = sph_harm(l_arr, m_arr, theta[:, None], phi[:, None])
    return 4 * np.pi * Y


def reconstruct_at_freq(freq_idx, freqs, positions, r_freq, N_cov, l_max=L_MAX_DEFAULT, n_modes=None):
    """Reconstruct the sky potential from one frequency bin's pulsar data.

    Parameters
    ----------
    freq_idx : int
        Index into `freqs` / the (2k, 2k+1) column pair of `r_freq`, `N_cov`.
    freqs, positions, r_freq, N_cov : as returned by `load_simulated_data`.
    l_max : int
        Maximum multipole included (l runs 2..l_max).
    n_modes : int or None
        SVD truncation rank; None keeps the full rank.

    Returns
    -------
    dict with keys: psi_lm, alm, dirty_map, fisher, singular_values, freq,
    plus fisher_cond_full (condition number of the untruncated Fisher matrix).
    """
    freq = freqs[freq_idx]
    n_psr = positions.shape[0]

    r_cos = r_freq[:, 2 * freq_idx + 1]
    r_sin = r_freq[:, 2 * freq_idx]
    r_data = r_cos - 1j * r_sin

    var_cos = N_cov[:, 2 * freq_idx + 1, 2 * freq_idx + 1]
    var_sin = N_cov[:, 2 * freq_idx, 2 * freq_idx]
    sigma2 = var_cos + var_sin
    inv_Cov = np.diag(1.0 / sigma2)

    lm_index = build_lm_index(l_max)
    Nl_vec = calc_Nl_vec(lm_index)
    theta, phi = hp.vec2ang(positions)

    A_f = 1.0 / (1j * 2 * np.pi * freq)
    D_raw = A_f * build_raw_design_matrix(theta, phi, lm_index)
    D = D_raw * (Nl_vec / 2)

    U, s_full, Vh = np.linalg.svd(D, full_matrices=False)

    Fisher_full = D.conj().T @ inv_Cov @ D
    fisher_cond_full = np.linalg.cond(Fisher_full)

    k = len(s_full) if n_modes is None else min(n_modes, len(s_full))
    U_k = U[:, :k]
    s_k = s_full[:k]
    V = Vh[:k, :].conj().T

    Delta = U_k @ np.diag(s_k)
    dirty_map = Delta.conj().T @ inv_Cov @ r_data
    Fisher = Delta.conj().T @ inv_Cov @ Delta
    Pi = np.linalg.inv(Fisher) @ dirty_map

    alm = V @ Pi
    psi_lm = alm * (Nl_vec / 2)

    return {
        "psi_lm": psi_lm,
        "alm": alm,
        "dirty_map": dirty_map,
        "fisher": Fisher,
        "singular_values": s_k,
        "freq": freq,
        "fisher_cond_full": fisher_cond_full,
        "n_psr": n_psr,
        "n_modes": k,
    }


def _print_summary(result):
    freq_nhz = result["freq"] * 1e9
    print(
        f"f = {freq_nhz:8.3f} nHz | "
        f"N_psr = {result['n_psr']:4d} | "
        f"N_basis = {len(result['psi_lm']):4d} | "
        f"n_modes = {result['n_modes']:4d} | "
        f"||P|| = {np.linalg.norm(result['psi_lm']):.4e} | "
        f"cond(Fisher_full) = {result['fisher_cond_full']:.4e}"
    )


def main(
    sample: str = typer.Option("psr_test_1", help="'psr_test_1' or 'psr_test_2'"),
    n_freqs: int = typer.Option(30, help="number of global Fourier modes"),
    l_max: int = typer.Option(L_MAX_DEFAULT, help="maximum multipole (l runs 2..l_max)"),
    n_modes: Optional[int] = typer.Option(None, help="SVD truncation rank; omit for full rank"),
    freq: Optional[int] = typer.Option(None, help="single frequency index; omit to run all"),
):
    if sample not in ("psr_test_1", "psr_test_2"):
        raise typer.BadParameter("sample must be 'psr_test_1' or 'psr_test_2'")

    data = load_simulated_data(sample, n_freqs=n_freqs)
    freqs, positions, r_freq, N_cov = data["freqs"], data["positions"], data["r_freq"], data["N_cov"]

    out_dir = Path("results")
    out_dir.mkdir(exist_ok=True)

    freq_indices = range(len(freqs)) if freq is None else [freq]
    results = []
    for k in freq_indices:
        result = reconstruct_at_freq(k, freqs, positions, r_freq, N_cov, l_max=l_max, n_modes=n_modes)
        _print_summary(result)
        np.savez(out_dir / f"sky_map_freq{k}.npz", **result)
        results.append(result)

    if freq is None:
        np.savez(
            out_dir / "sky_map_all.npz",
            freqs=freqs,
            psi_lm=np.array([r["psi_lm"] for r in results]),
            alm=np.array([r["alm"] for r in results]),
            dirty_map=np.array([r["dirty_map"] for r in results]),
            fisher=np.array([r["fisher"] for r in results]),
            singular_values=np.array([r["singular_values"] for r in results]),
            fisher_cond_full=np.array([r["fisher_cond_full"] for r in results]),
        )


if __name__ == "__main__":
    typer.run(main)
