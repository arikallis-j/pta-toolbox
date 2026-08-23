import jax
jax.config.update("jax_enable_x64", True)

import numpy as np
import s2fft

def mat2vec_sph(F_lm, nside):
    L = 2 * nside
    f_lm = s2fft.sampling.s2_samples.flm_2d_to_1d(F_lm, L)
    return f_lm

def vec2mat_sph(f_lm, nside):
    L = 2 * nside
    F_lm = s2fft.sampling.s2_samples.flm_1d_to_2d(f_lm, L)
    return F_lm

def calc_Nl(nside):
    L = 2 * nside
    nl = np.zeros(L**2, dtype=np.float64)
    idx = 0 
    for l in range(L):
        n_m = 2 * l + 1
        if l>1:
            nl[idx : idx + n_m] = np.sqrt(2/((l-1)*(l)*(l+1)*(l+2)))
        else:
            nl[idx : idx + n_m] = np.nan
        idx += n_m
    Nl = vec2mat_sph(nl, nside)
    return Nl

def lr2gc_alm(alm_l, alm_r):
    alm_g = alm_l + alm_r
    alm_c = 1/(1j) * (alm_l - alm_r)
    return alm_g, alm_c

def gc2lr_alm(alm_g, alm_c):
    alm_l = 1/2 * (alm_g + 1j * alm_c)
    alm_r = 1/2 * (alm_g - 1j * alm_c)
    return alm_l, alm_r

def alm2plm(alm_g, alm_c, Nl):
    plm_g = alm_g * Nl / 2
    plm_c = alm_c * Nl / 2
    return plm_g, plm_c

def plm2alm(plm_g, plm_c, Nl):
    alm_g = 2 / Nl * plm_g 
    alm_c = 2 / Nl * plm_c
    return alm_g, alm_c

def generate_flm(nside, key=42, L_lower=0):
    L = 2 * nside
    rng = np.random.default_rng(key)
    F_lm = 1/np.sqrt(2) * s2fft.utils.signal_generator.generate_flm(rng, L, L_lower=L_lower)
    return F_lm

def generate_complex_normal(nside, var=1.0, key=42):
    Npix = 12 * nside**2
    rng = np.random.default_rng(key)
    f = s2fft.utils.signal_generator.complex_normal(rng, Npix, var)
    return f

def forward_sfft2(f, nside, niter=15, spin=0, L_lower=None):
    L = 2 * nside
    if L_lower is None:
        L_lower = abs(spin)
    F_lm = s2fft.forward(f, L, spin=spin, L_lower=L_lower, nside=nside, iter=niter, sampling="healpix", method="numpy")
    return F_lm

def inverse_sfft2(F_lm, nside, spin=0, L_lower=None):
    L = 2 * nside
    if L_lower is None:
        L_lower = abs(spin)
    f = s2fft.inverse(F_lm, L, spin=spin, L_lower=L_lower, nside=nside, sampling="healpix", method="numpy")
    return f