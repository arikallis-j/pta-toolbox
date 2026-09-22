import ptatoolbox as pta
import xarray as xr
import typer

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from scipy.special import sph_harm_y
from matplotlib.colors import LogNorm

def build_lm_index(l_max: int):
    """List of (l, m) pairs for l = 0..l_max, m = -l..l."""
    return [(l, m) for l in range(0, l_max + 1) for m in range(-l, l + 1)]

def calc_Nl_vec(lm_index):
    """Nl(l) for l >= 2, zero for l < 2."""
    l_arr = np.array([l for l, _ in lm_index], dtype=np.float64)
    Nl = np.where(
        l_arr >= 2,
        np.sqrt(2.0 / ((l_arr - 1) * l_arr * (l_arr + 1) * (l_arr + 2))),
        0.0
    )
    return Nl

def truncate_modes(s, Fisher, rtol=1e-2, return_mask=True):
    """
    Combined truncation: singular value criterion + Fisher condition number.
    
    Parameters
    ----------
    s : ndarray, shape (n_psr,)
        Singular values from SVD of D.
    Fisher : ndarray, shape (n_freqs, n_psr, n_psr)
        Fisher matrix per frequency.
    rtol : float
        Relative threshold for singular values.
    
    Returns
    -------
    mask : ndarray, shape (n_psr,), bool
        True for modes to keep.
    """
    # геометрический критерий (не зависит от частоты)
    mask_svd = s / s.max() > rtol
    
    # критерий по числу обусловленности (усредняем по частотам)
    # cond = np.array([
    #     np.linalg.cond(Fisher[k]) for k in range(len(Fisher))
    # ])
    # print(f"Fisher condition numbers: min={cond.min():.2e}, max={cond.max():.2e}")
    return mask_svd

def build_raw_design_matrix(theta, phi, lm_index):
    """D_raw[i, nu] = 4*pi * Y_lm(theta_i, phi_i). Shape: (n_psr, n_modes)."""
    l_arr = np.array([l for l, _ in lm_index])[None, :]
    m_arr = np.array([m for _, m in lm_index])[None, :]
    Y = np.where(
        l_arr >= 2,
        sph_harm_y(l_arr, m_arr, theta[:, None], phi[:, None]),
        0.0
    )
    return 4 * np.pi * Y

def make_sky(ds_obs, lmax, name='sample'):
    lm_index = build_lm_index(lmax)
    Nl_vec = calc_Nl_vec(lm_index)

    theta, phi, freq = ds_obs.theta.values, ds_obs.phi.values, ds_obs.freq.values
    r_freq, N_cov = ds_obs.r_freq.values, ds_obs.N_cov.values

    inv_Cov = np.linalg.inv(N_cov)
    A_f = 1.0 / (1j * 2 * np.pi * freq)
    D_raw = build_raw_design_matrix(theta, phi, lm_index)
    D = D_raw * (Nl_vec / 2)

    U, s, Vh = np.linalg.svd(D, full_matrices=False)
    Sigma = np.diag(s)
    Delta = U @ Sigma
    V = Vh.T.conj()
    Delta_f = A_f[:, None, None] * Delta[None, :, :]  # shape (n_freqs, n_psr, n_psr)

    # Fisher matrix: (n_freqs, n_psr, n_psr)
    Fisher = np.einsum('kia,kib->kab', Delta_f.conj(), np.einsum('kij,kjb->kib', inv_Cov, Delta_f))
    # dirty map: (n_freqs, n_psr)
    X = np.einsum('kia,ki->ka', Delta_f.conj(), np.einsum('kij,kj->ki', inv_Cov, r_freq))
    # clean map: (n_freqs, n_psr)
    P = np.einsum('kab,kb->ka', np.linalg.inv(Fisher) , X)  # (n_freqs, n_psr)

    l_arr = np.array([l for l, _ in lm_index])
    m_arr = np.array([m for _, m in lm_index])
    mode_idx = np.arange(len(s))

    ds_sky = xr.Dataset(
        {
            "P":  (["freq", "mode"], P),              # complex128
            "F": (["freq", "mode_i", "mode_j"], Fisher),  # complex128
            "V":   (["harm", "mode"], V),             # complex128
            "s":   (["mode"], s),                     # float64
        },
        coords={
            "freq": freq,
            "mode": mode_idx,
            "mode_i": mode_idx,
            "mode_j": mode_idx,
            "harm": np.arange(len(lm_index)),
            "l":    (["harm"], l_arr),
            "m":    (["harm"], m_arr),
        },
        attrs={
            "name": name,
            "lmax": lmax,
        }
    )
    return ds_sky

def plot_fisher(ds_sky, freq_idx=0):
    Fisher_k = ds_sky.F.values[freq_idx]  # (n_psr, n_psr)
    freq_nHz = ds_sky.freq.values[freq_idx] / 1e-9

    fig, ax = plt.subplots()
    im = ax.pcolormesh(
        np.abs(Fisher_k),
        cmap="viridis",
        norm=LogNorm()
    )
    plt.colorbar(im, ax=ax)
    ax.set_title(f"Fisher matrix |F| at f = {freq_nHz:.2f} nHz")
    ax.set_xlabel("mode j")
    ax.set_ylabel("mode i")
    plt.show()

def plot_IV(I, V, name='stokes-I', show=True, path=None, grid=True, vmin=0.1):
    plt.figure(figsize=(6, 8))
    I0, V0 = np.max(I), np.max(np.abs(V))
    I_plot = np.clip(I / I0, vmin, None)
    
    hp.mollview(I_plot, sub=211, title='Total Intensity',
                cmap='plasma', norm='log', min=vmin, max=1.0)

    if grid:
        hp.graticule()
    hp.mollview(V, sub=212, title='Total Сhirality', cmap='coolwarm', min=-V0, max=V0)
    if grid:
        hp.graticule()
    plt.tight_layout()
    if path is not None:
        plt.savefig(path / f"{name}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def main(nside: int = 50, n_freqs: int = -1, rtol: float = 1e-2, sample: str = "psr_test_1", show: bool = False):
    lmax = 2 * nside
    dm = pta.utils.DataManager()
    path = dm.create_experiment(f"obs-{sample}")
    datapath = dm.raw / sample
    psi_path = dm.make_dir(path, "psi_maps")
    h_path = dm.make_dir(path, "h_maps")
    I_path = dm.make_dir(path, "I_maps")

    # pulsars = pta.utils.load_pulsars(path=datapath)
    
    # ds_obs = pta.utils.make_observations(pulsars, name=sample, n_freqs=n_freqs)
    # ds_obs.to_netcdf(path / f"{sample}-obs.nc", engine="h5netcdf")
    with xr.open_dataset(path / f"{sample}-obs.nc" , engine="h5netcdf") as data:
        ds_obs = data.load()

    # ds_sky = make_sky(ds_obs, lmax, name=sample)
    # ds_sky.to_netcdf(path / f"{sample}-sky.nc", engine="h5netcdf")
    with xr.open_dataset(path / f"{sample}-sky.nc" , engine="h5netcdf") as data:
        ds_sky = data.load()

    print(ds_sky)


    lm_index = build_lm_index(lmax)
    Nl_vec = calc_Nl_vec(lm_index)

    mask = truncate_modes(ds_sky.s.values, ds_sky.F.values, rtol=rtol)
    P = ds_sky.P.values[:, mask]       # (n_freqs, n_modes_kept)
    V  = ds_sky.V.values[:, mask] 
    print(f"P.shape {P.shape}") 

    alm = np.einsum('na,ka->kn', V, P)
    print(f"alm.shape {alm.shape}")

    plm = alm * (Nl_vec / 2)
    print(f"plm.shape {plm.shape}")

    # psi_maps = []

    # for k in range(len(ds_sky.freq)):
    #     Plm = pta.vec2mat_sph(plm[k], nside)
    #     h_plus, h_cross = pta.psi2hpol(nside, Plm)
    #     pta.plot_intensities(h_plus, h_cross, name=f'I_{ds_sky.freq[k]/1e-9:.2f}nHz', show=show, path=I_path, grid=True, log=True)



    delta_f = ds_sky.freq.values[1] - ds_sky.freq.values[0]  # 1/tspan, равномерная сетка
    print(f"delta_f = {delta_f/1e-9:.2f} nHz")
    
    I_total, V_total, Plm_total = None, None, None
    for k in range(len(ds_sky.freq)):
        Plm = pta.vec2mat_sph(plm[k], nside)
        h_plus, h_cross = pta.psi2hpol(nside, Plm)
        
        I_k = np.abs(h_plus)**2 + np.abs(h_cross)**2
        I_k, _, _, V_k = pta.Stokes_parameters(h_plus, h_cross)
        if I_total is None:
            I_total = np.zeros_like(I_k)
        if V_total is None:
            V_total = np.zeros_like(V_k)
        I_total += I_k * delta_f
        V_total += V_k * delta_f

        # if Plm_total is None:
        #     Plm_total = np.zeros_like(Plm)
        # Plm_total += Plm * delta_f

        plot_IV(I_k, V_k, name=f'I_{ds_sky.freq[k]/1e-9:.2f}nHz', show=show, path=I_path, grid=True)
    
    plot_IV(I_total, V_total, name=f'IV_total', show=show, path=path, grid=True)

if __name__ == '__main__':
    typer.run(main)

    # for k in range(2):#range(len(ds_sky.freq)):
    #     Plm = pta.vec2mat_sph(plm[k], nside)
    #     psi = pta.inverse_sfft2(Plm, nside, spin=0)
    #     pta.plot_map(psi, name=f'psi_{ds_sky.freq[k]/1e-9:.2f}nHz', path=psi_path, show=show)
    # for k in range(2):#range(len(ds_sky.freq)):
    #     Plm = pta.vec2mat_sph(plm[k], nside)
    #     h_plus, h_cross = pta.psi2hpol(nside, Plm)
    #     pta.plot_map(h_plus, name=f'h+_{ds_sky.freq[k]/1e-9:.2f}nHz', path=h_path, show=show)
    #     pta.plot_map(h_cross, name=f'hx_{ds_sky.freq[k]/1e-9:.2f}nHz', path=h_path, show=show)
    