import ptatoolbox as pta
import numpy as np
import healpy as hp
import typer
import s2fft

import matplotlib.pyplot as plt

def make_delta_map(nside, ra_s, dec_s, normalize=False):
    # Проверка корректности nside
    if not hp.isnsideok(nside):
        raise ValueError("nside должно быть степенью двойки (например, 1, 2, 4, 8, ...)")

    npix = hp.nside2npix(nside)
    m = np.zeros(npix, dtype=np.float64)

    theta = np.radians(90.0 - dec_s)
    phi = np.radians(ra_s)

    ipix = hp.ang2pix(nside, theta, phi)

    if normalize:
        pix_area = hp.nside2pixarea(nside)
        m[ipix] = 1.0 / pix_area
    else:
        m[ipix] = 1.0
    return m


def main(nside: int = 50, real: bool = False, lobs: int = 50, polar: str = '+'):
    lmax = 2*nside
    dm = pta.utils.DataManager()
    path = dm.create_experiment("point")

    ra, dec, A = 0.0, 0.0, 1.0
    Nl = pta.calc_Nl(nside)
    delta = make_delta_map(nside, ra, dec, normalize=True)
    
    delta_lm_left = pta.forward_sfft2(delta, nside, spin=+2)
    delta_lm_right = pta.forward_sfft2(delta, nside, spin=-2)
    delta_left = pta.inverse_sfft2(delta_lm_left, nside, spin=+2)
    delta_right = pta.inverse_sfft2(delta_lm_right, nside, spin=-2)

    if polar == '+':
        p_left, p_right = A * delta_left, A * delta_right
        p_plus, p_cross = pta.circ2lin(p_left, p_right)
        p_plus, p_cross = p_plus, np.zeros_like(p_cross)
        p_left, p_right = pta.lin2circ(p_plus, p_cross)
    elif polar == 'x':
        p_left, p_right = A * delta_left, - A * delta_right
        p_plus, p_cross = pta.circ2lin(p_left, p_right)
        p_plus, p_cross = np.zeros_like(p_plus), p_cross
        p_left, p_right = pta.lin2circ(p_plus, p_cross)
    elif polar == 'l':
        p_left, p_right = 2*A * delta_left, np.zeros_like(delta_right)
        p_plus, p_cross = pta.circ2lin(p_left, p_right)
    elif polar == 'r':
        p_left, p_right = np.zeros_like(delta_left), 2*A * delta_right
        p_plus, p_cross = pta.circ2lin(p_left, p_right)
    else:
        pass
    
    pta.plot_map(p_plus, name='h-plus', path=path, show=False)
    pta.plot_map(p_cross, name='h-cross', path=path, show=False)
    pta.plot_intensities(p_plus, p_cross, path=path, show=False, grid=False)
    pta.plot_polarization_degrees(p_plus, p_cross, path=path, show=False)
    pta.plot_Stokes_parameters(p_plus, p_cross, path=path, show=False)

    alm_l = pta.forward_sfft2(p_left, nside, spin=+2)
    alm_r = pta.forward_sfft2(p_right, nside, spin=-2)
    
    alm_g, alm_c = pta.lr2gc_alm(alm_l, alm_r)
    plm_g, plm_c = pta.alm2plm(alm_g, alm_c, Nl)

    psi_g = pta.inverse_sfft2(plm_g, nside, spin=0)
    psi_c = pta.inverse_sfft2(plm_c, nside, spin=0)
    
    pta.plot_map(psi_g, name='psi-g_aph', path=path, show=False)
    pta.plot_map(psi_c, name='psi-c_aph',path=path, show=False)
    pta.plot_maps(psi_g, name='psi-g_ri', path=path, show=False)
    pta.plot_maps(psi_c, name='psi-c_ri',path=path, show=False)

if __name__ == '__main__':
    typer.run(main)