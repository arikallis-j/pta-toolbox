import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import ptatoolbox as pta
import typer

def main(exp: str = 'test', nside: int = 20, flat: bool = False, proj: str = 'ri'):
    dm = pta.utils.DataManager()
    sim_path = dm.make_dir(path=dm.raw, name=exp)
    res_path = dm.make_dir(path=dm.processed, name=exp)

    lmax = 2 * nside
    npix = hp.nside2npix(nside)
    dOmega = 4 * np.pi / npix

    a_lm = np.ones(lmax**2, dtype=complex)
    A_lm = pta.vec2mat(a_lm, lmax)
    a_lm_est = pta.mat2vec(A_lm, lmax)
    print(np.max(np.abs(a_lm - a_lm_est)))

    a_lm = pta.generate_flm(lmax, key=42, lmin=2)
    nl = pta.generate_nl(lmax)
    ell = pta.generate_ell(lmax)
    p_lm = nl * a_lm
    a_lm_est = np.where(ell > 1, p_lm/nl, 0.0)

    map_pix = pta.harm2pix(p_lm, lmax)
    
    if flat:
        fig, axis = plt.subplots(2, figsize=(10, 10))
    else:
        fig, axis = plt.subplots(2, subplot_kw={'projection': 'mollweide'}, figsize=(10, 10))

    if proj == 'ri':
        # real/imag maps
        sky_real = pta.plot_sky_real(map_pix, axis[0], cmap='coolwarm', hours=True)
        sky_imag = pta.plot_sky_imag(map_pix, axis[1], cmap='coolwarm', hours=True)
        pta.plot_colorbar_healpix(sky_real, fig, axis[0], orientation='vertical', shrink=0.9, aspect=20, pad=0.05)
        pta.plot_colorbar_healpix(sky_imag, fig, axis[1], orientation='vertical', shrink=0.9, aspect=20, pad=0.05)
    elif proj == 'ap':
        # amp/phase maps
        sky_real = pta.plot_sky_amplitude(map_pix, axis[0], cmap='viridis', hours=True)
        sky_imag = pta.plot_sky_phase(map_pix, axis[1], hours=True)
        pta.plot_colorbar_healpix(sky_real, fig, axis[0], orientation='vertical', shrink=0.9, aspect=20, pad=0.05)
        pta.plot_colorbar_healpix(sky_imag, fig, axis[1], orientation='vertical', shrink=0.9, aspect=20, pad=0.05)

    fig.suptitle("Gradient Potential $\Psi^G$")
    fig.tight_layout()
    plt.show()

    return None

if __name__ == '__main__':
    typer.run(main)