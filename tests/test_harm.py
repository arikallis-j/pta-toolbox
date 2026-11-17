import matplotlib.pyplot as plt
import numpy as np
import ptatoolbox as pta
import typer

def main(exp: str = 'test-harm', nside: int = 20):
    dm = pta.utils.DataManager()
    sim_path = dm.make_dir(path=dm.raw, name=exp)
    res_path = dm.make_dir(path=dm.processed, name=exp)

    lmax = 2 * nside

    a_lm = np.ones(lmax**2, dtype=complex)
    A_lm = pta.vec2mat(a_lm, lmax)
    a_lm_est = pta.mat2vec(A_lm, lmax)
    print(np.max(np.abs(a_lm - a_lm_est)))

    a_lm = pta.generate_flm(lmax, key=42, lmin=2)
    a_pix = pta.harm2pix(a_lm, lmax)
    a_lm_est = pta.pix2harm(a_pix, lmax)
    print(np.max(np.abs(a_lm - a_lm_est)))
    
    nl = pta.generate_nl(lmax)
    ell = pta.generate_ell(lmax)
    p_lm = nl * a_lm
    a_lm_est = np.where(ell > 1, p_lm/nl, 0.0)
    print(np.max(np.abs(a_lm - a_lm_est)))

    Cl = 0.0
    n = 1000
    for k in range(n):
        a_lm = pta.generate_flm(lmax, key=k, lmin=0)
        cl = pta.angular_spectrum(a_lm, lmax)
        Cl += cl
    cl = Cl/n
    print(np.max(np.abs(cl - 1.0)))

    return None

if __name__ == '__main__':
    typer.run(main)