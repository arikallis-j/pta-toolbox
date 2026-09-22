import ptatoolbox as pta
import pta_catalog as ptacat
import numpy as np
import typer

from scipy.special import sph_harm_y as sph_harm


def plot_eigenmodes(V, s, nside, Nl_mat, path=None, show=False):
    """Print singular values and plot the spatial (psi) pattern of each right-singular eigenmode of D."""
    print(f"Singular values: {s}")
    for k in range(V.shape[1]):
        mode_alm = pta.vec2mat_sph(V[:, k], nside)
        mode_plm = Nl_mat / 2 * mode_alm
        mode_psi = pta.inverse_sfft2(mode_plm, nside, spin=0)
        pta.plot_maps(mode_psi, name=f'mode-{k+1}', path=path, show=show)

def build_design_matrix(theta, phi, lmax):
    """Raw (unnormalized) design matrix: D[i, idx] = 4*pi*Y_lm(theta_i, phi_i), zero for l < 2."""
    n_psr = len(theta)
    D = np.zeros((n_psr, lmax**2), dtype=np.complex128)
    for i in range(n_psr):
        th, ph = theta[i], phi[i]
        idx = 0
        for l in range(lmax):
            for m in range(-l, l + 1):
                if l >= 2:
                    D[i, idx] = 4 * np.pi * sph_harm(l, m, th, ph)
                idx += 1
    return D

SPHERE = {
    'coords': 'sphere',
    'params': {
        'seed_coord': 42,
        'radius': 1.0,
    }
}

CAP = {
    'coords': 'cap',
    'params': {
        'seed_coord': 42,
        'radius': 1.0,
        'ra_0': 0.0, 
        'dec_0': 0.0, 
        'alpha': 30.0,
    }
}

def main(nside: int = 50, n_psr: int = 10, n_modes: int = None, sigma: float = 0.1, show: bool = False):
    lmax = 2 * nside

    dm = pta.utils.DataManager()
    path = dm.create_experiment("research")
    plm_g = pta.make_stochastic(nside, key=42)
    plm_vec = pta.mat2vec_sph(plm_g, nside)
    psrcat = ptacat.make_catalog(n_psr, config=SPHERE)
    ptacat.plot_catalog(psrcat, path=path, show=show)

    data = psrcat.data
    ra_rad = np.deg2rad(data['RAJD'].values)
    dec_rad = np.deg2rad(data['DECJD'].values)
    theta = np.pi / 2 - dec_rad
    phi = ra_rad

    D_raw = build_design_matrix(theta, phi, lmax)

    r_true = D_raw @ plm_vec
    noise = 1/np.sqrt(2) * (np.random.normal(0, sigma, n_psr) + 1j * np.random.normal(0, sigma, n_psr))
    r_data = r_true + noise

    Nl_mat = pta.calc_Nl(nside)
    Nl_vec = pta.mat2vec_sph(Nl_mat, nside)
    D = D_raw * (Nl_vec / 2)

    Cov = sigma**2 * np.identity(n_psr)
    inv_Cov = np.linalg.inv(Cov)

    U, s, Vh = np.linalg.svd(D, full_matrices=False)

    if n_modes is not None:
        k = n_modes
    else:
        k = len(s)
    print(f"Ранг: {k} из {len(s)}")

    U = U[:, :k]
    s = s[:k]
    Vh = Vh[:k, :]

    Sigma = np.diag(s)
    Delta = U @ Sigma
    V = Vh.T.conj()

    if n_modes is not None:
        plot_eigenmodes(V, s, nside, Nl_mat, path=path, show=show)

    Chi = Delta.T.conj() @ inv_Cov @ r_data
    Fisher = Delta.T.conj() @ inv_Cov @ Delta
    Pi = np.linalg.inv(Fisher) @ Chi
    a_vec_est = V @ Pi
    alm_est = pta.vec2mat_sph(a_vec_est, nside)
    plm_est = alm_est * Nl_mat/2
    psi_est = pta.inverse_sfft2(plm_est, nside, spin=0)
    psi_g = pta.inverse_sfft2(plm_g, nside, spin=0)

    pta.plot_map(psi_g, name='psi_g', path=path, show=show)
    pta.plot_map(psi_est, name=f'psi_est', path=path, show=show)


if __name__ == '__main__':
    typer.run(main)