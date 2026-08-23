import ptatoolbox as pta
import numpy as np
import healpy as hp
import typer
import s2fft

import matplotlib.pyplot as plt

def main(nside: int = 50, real: bool = False, lobs: int = 50):
    lmax = 2*nside
    dm = pta.utils.DataManager()
    path = dm.create_experiment("stochastic")

    Nl = pta.calc_Nl(nside)
    plm_g = Nl/2 * pta.generate_flm(nside, key=42)
    plm_c = Nl/2 * pta.generate_flm(nside, key=43)

    if real:
        plm_g, plm_c = plm_g, np.zeros_like(plm_c)

    if lobs != 0:
        plm_g_vec = pta.mat2vec_sph(plm_g, nside)
        plm_c_vec = pta.mat2vec_sph(plm_c, nside)
        plm_g_vec[lobs**2::] = 0.0
        plm_c_vec[lobs**2::] = 0.0
        plm_g = pta.vec2mat_sph(plm_g_vec, nside)
        plm_c = pta.vec2mat_sph(plm_c_vec, nside)

    psi_g = pta.inverse_sfft2(plm_g, nside, spin=0)
    psi_c = pta.inverse_sfft2(plm_c, nside, spin=0)
    pta.plot_map(psi_g, name='psi-g', path=path, show=False)
    pta.plot_map(psi_c, name='psi-c',path=path, show=False)

    alm_g, alm_c = pta.plm2alm(plm_g, plm_c, Nl)
    alm_l, alm_r = pta.gc2lr_alm(alm_g, alm_c)

    h_left = pta.inverse_sfft2(alm_l, nside, spin=2)
    h_right = pta.inverse_sfft2(alm_r, nside, spin=-2)
    h_plus, h_cross = pta.circ2lin(h_left, h_right) 

    pta.plot_map(h_plus, name='h-plus', path=path, show=False)
    pta.plot_map(h_cross, name='h-cross', path=path, show=False)
    pta.plot_intensities(h_plus, h_cross, path=path, show=False)

if __name__ == '__main__':
    typer.run(main)