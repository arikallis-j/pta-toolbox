import jax
jax.config.update("jax_enable_x64", True)

import numpy as np
import healpy as hp
import s2fft

def forward_sfft2(f, lmax, spin=0, L_lower=None, niter=15):
    if L_lower is None:
        L_lower = abs(spin)
    F_lm = s2fft.forward(f, lmax, spin=spin, L_lower=L_lower, nside=lmax//2, iter=niter, sampling="healpix", method="numpy")
    return F_lm

def inverse_sfft2(F_lm, lmax, spin=0, L_lower=None):
    if L_lower is None:
        L_lower = abs(spin)
    f = s2fft.inverse(F_lm, lmax, spin=spin, L_lower=L_lower, nside=lmax//2, sampling="healpix", method="numpy")
    return f

def generate_nl(lmax):
    nl = np.zeros(lmax**2, dtype=np.float64)
    idx = 0 
    for l in range(lmax):
        n_m = 2 * l + 1
        if l>1:
            nl[idx : idx + n_m] = np.sqrt(2/((l-1)*(l)*(l+1)*(l+2)))
        else:
            nl[idx : idx + n_m] = 0.0
        idx += n_m
    return nl

def generate_ell(lmax):
    ell = np.zeros(lmax**2, dtype=np.float64)
    idx = 0 
    for l in range(lmax):
        n_m = 2 * l + 1
        ell[idx : idx + n_m] = l
        idx += n_m
    return ell

def generate_flm(lmax, key=42, lmin=0):
    rng = np.random.default_rng(key)
    ell = generate_ell(lmax)
    f_lm = 1/np.sqrt(2) * (rng.normal(0.0, 1.0, size=lmax**2) + 1j * rng.normal(0.0, 1.0, lmax**2))
    f_lm[ell<lmin] = 0.0 
    return f_lm

def angular_spectrum(f_lm, lmax):
    cl = np.zeros(lmax, dtype=np.float64)
    idx = 0
    for l in range(lmax):
        n_m = 2 * l + 1
        flm_l = f_lm[idx : idx + n_m]
        cl[l] = np.mean(np.abs(flm_l) ** 2)
        idx += n_m
    return cl

def mat2vec(F_lm, lmax):
    f_lm = s2fft.sampling.s2_samples.flm_2d_to_1d(F_lm, lmax)
    return f_lm

def vec2mat(f_lm, lmax):
    F_lm = s2fft.sampling.s2_samples.flm_1d_to_2d(f_lm, lmax)
    return F_lm

def pix2harm(f_pix, lmax, spin=0):
    F_lm = forward_sfft2(f_pix, lmax, spin=spin)
    f_lm = mat2vec(F_lm, lmax)
    return f_lm

def harm2pix(f_lm, lmax, spin=0):
    F_lm = vec2mat(f_lm, lmax)
    f_pix = inverse_sfft2(F_lm, lmax, spin=spin)
    return f_pix