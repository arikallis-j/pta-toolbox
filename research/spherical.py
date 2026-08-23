import ptatoolbox as pta
import typer

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from time import perf_counter

def main(nside: int = 100, niter: int = 15, ndeg: int = 10, spin: int = 0, plot: bool = False):
    F_lm_data = pta.generate_flm(nside, L_lower=spin)
    f_data = pta.inverse_sfft2(F_lm_data, nside, spin=spin)
    
    start = perf_counter()
    F_lm = pta.forward_sfft2(f_data, nside, spin=spin, niter=niter)
    end = perf_counter()
    print(f"forward ~ {end - start:.2f} sec")

    start = perf_counter()
    f_alpha = pta.inverse_sfft2(F_lm, nside, spin=spin)
    end = perf_counter()
    print(f"inverse ~ {end - start:.2f} sec")

    pta.plot_map(f_data)
    pta.plot_map(f_alpha)

    start = perf_counter()
    f_lm = pta.mat2vec_sph(F_lm, nside)
    end = perf_counter()
    print(f"mat2vec ~ {end - start:.2f} sec")
    print(f"f_lm ~ {f_lm.shape}")

    start = perf_counter()
    F_lm = pta.vec2mat_sph(f_lm, nside)
    end = perf_counter()
    print(f"vec2mat ~ {end - start:.2f} sec")
    print(f"F_lm ~ {F_lm.shape}")

    error = np.abs(f_data - f_alpha)/np.abs(f_data)
    print(f"relative error ~ {error.max():.2e}")
    
    if plot:
        hp.mollview(error, cmap="viridis")
        plt.show()

    for k in range(ndeg):
        F_lm = pta.forward_sfft2(f_alpha, nside, spin=spin, niter=niter)
        f_alpha = pta.inverse_sfft2(F_lm, nside, spin=spin)
        error = np.abs(f_data - f_alpha)/np.abs(f_data)
        print(f"relative error deg[{k+1}] ~ {error.max():.2e}")

if __name__ == '__main__':
    typer.run(main)