import ptatoolbox as pta
import numpy as np
import healpy as hp
import typer
import s2fft

import matplotlib.pyplot as plt

def main(nside: int = 50):
    lmax = 2*nside
    dm = pta.utils.DataManager()
    path = dm.create_experiment("modes")

    Nl = pta.calc_Nl(nside)
    alm_left = pta.generate_flm(nside, key=42, L_lower=2)
    alm_right = pta.generate_flm(nside, key=43, L_lower=2)

    h_left = pta.inverse_sfft2(alm_left, nside, spin=2)
    h_right = pta.inverse_sfft2(alm_right, nside, spin=-2)
    h_plus, h_cross = pta.circ2lin(h_left, h_right) 

    alm_g, alm_c = pta.lr2gc_alm(alm_left, alm_right)
    plm_g, plm_c = pta.alm2plm(alm_g, alm_c, Nl)

    psi_g = pta.inverse_sfft2(plm_g, nside, spin=0)
    psi_c = pta.inverse_sfft2(plm_c, nside, spin=0)

    cl = pta.angular_spectrum(plm_g, nside)
    pta.plot_angular_spectrum(cl, lmax, log=True)

    pta.plot_map(psi_g)
    pta.plot_map(psi_c)

if __name__ == '__main__':
    typer.run(main)