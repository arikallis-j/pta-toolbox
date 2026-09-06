import ptatoolbox as pta
import pta_catalog as ptacat
import numpy as np
import healpy as hp
import typer
import s2fft

import matplotlib.pyplot as plt
from scipy.special import sph_harm_y as sph_harm


def main(nside: int = 50, n_psr: int = 10, show: bool = False, sigma: float = 0.1, modes: int = None, lobs: int = None):
    lmax = 2 * nside
    dm = pta.utils.DataManager()
    path = dm.create_experiment("eigenmodes")

    # --- Генерация исходных данных ---
    Nl = pta.calc_Nl(nside)
    Plm_true = Nl / 2 * pta.generate_flm(nside, key=42)
    P_true = pta.inverse_sfft2(Plm_true, nside, spin=0)
    pta.plot_map(P_true, name='P_true', path=path, show=show)

    # alm_g, alm_c = pta.plm2alm(Plm_true, np.zeros_like(Plm_true), Nl)
    # alm_l, alm_r = pta.gc2lr_alm(alm_g, alm_c)

    # h_left = pta.inverse_sfft2(alm_l, nside, spin=2)
    # h_right = pta.inverse_sfft2(alm_r, nside, spin=-2)
    # h_plus, h_cross = pta.circ2lin(h_left, h_right) 

    # pta.plot_map(h_plus, name='h-plus-true', path=path, show=False)
    # pta.plot_map(h_cross, name='h-cross-true', path=path, show=False)
    # pta.plot_intensities(h_plus, h_cross, name='stokes-I-true', path=path, show=False)


    # Преобразуем в вектор и обрезаем l=0,1
    plm_full = pta.mat2vec_sph(Plm_true, nside)
    plm_trimmed = plm_full[4:]   # убираем первые 4 моды (l=0,1)
    n_modes = lmax**2 - 4

    # --- Каталог пульсаров ---
    psrcat = ptacat.make_catalog(n_psr, config={
    # 'coords': 'cap',
    # 'params': {
    #     'seed_coord': 42,
    #     'ra_0': 0.0, 
    #     'dec_0': 0.0, 
    #     'radius': 30.0,
    #     }
    })
    ptacat.plot_catalog(psrcat, path, show=show)

    data = psrcat.data
    ra_rad = np.deg2rad(data['RAJD'].values)
    dec_rad = np.deg2rad(data['DECJD'].values)
    theta = np.pi / 2 - dec_rad   # полярный угол
    phi = ra_rad                  # азимутальный угол
    
    # --- Построение матрицы отклика D (с 4π) ---
    D = np.zeros((n_psr, n_modes), dtype=np.complex128)
    for i in range(n_psr):
        th, ph = theta[i], phi[i]
        idx = 0
        for l in range(2, lmax):          # только l >= 2
            for m in range(-l, l + 1):
                D[i, idx] = 4 * np.pi * sph_harm(l, m, th, ph)   # правильный порядок: (l, m, theta, phi)
                idx += 1

    # --- Данные (отклик) с шумом ---
    r_true = D @ plm_trimmed
    noise = np.random.normal(0, sigma, n_psr) + 1j * np.random.normal(0, sigma, n_psr)
    r_data = r_true + noise

    Cov = sigma**2 * np.identity(n_psr)
    inv_Cov = np.linalg.inv(Cov)

    # --- SVD и восстановление ---
    U, s, Vh = np.linalg.svd(D, full_matrices=False)
    # Отсекаем малые сингулярные числа (по порогу или по энергии)

    
    if modes is not None:
        k = modes
    else:
        k = len(s[s>1])
    print(f"Ранг: {k} из {len(s)}")

    U_k = U[:, :k]
    Sigma_k = np.diag(s[:k])
    Vh_k = Vh[:k, :]

    Delta = U_k @ Sigma_k
    V = Vh_k.T.conj()

    print(V.shape)
    alpha = 0
    eigenv_trimmed_lm = V[:, alpha]
    eigenv_full_lm = np.zeros(lmax**2, dtype=complex)
    eigenv_full_lm[4:] = eigenv_trimmed_lm 
    Eigenv_lm = pta.vec2mat_sph(eigenv_full_lm, nside)
    Eigenv = pta.inverse_sfft2(Eigenv_lm, nside, spin=0)
    pta.plot_map(Eigenv, name='Eigenv', path=path, show=show)


    Chi = Delta.T.conj() @ inv_Cov @ r_data
    Fisher = Delta.T.conj() @ inv_Cov @ Delta
    Pi = np.linalg.inv(Fisher) @ Chi

    plm_trimmed_est = V @ Pi

    plm_trimmed_obs = V @ Vh_k @ plm_trimmed 
    # # --- Восстановление полного вектора ---
    plm_full_obs = np.zeros(lmax**2, dtype=complex)
    plm_full_obs[4:] = plm_trimmed_obs 
    if lobs is not None:
        plm_full_obs[lobs**2+1::] = 0.0
    Plm_obs = pta.vec2mat_sph(plm_full_obs, nside)
    P_obs = pta.inverse_sfft2(Plm_obs, nside, spin=0)
    pta.plot_map(P_obs, name='P_obs', path=path, show=show)


    # # # --- Восстановление полного вектора ---
    # plm_full_est = np.zeros(lmax**2, dtype=complex)
    # plm_full_est[4:] = plm_trimmed_est 

    # # Преобразуем в карту и визуализируем
    # Plm_est = pta.vec2mat_sph(plm_full_est, nside)

    # alm_g, alm_c = pta.plm2alm(Plm_est, np.zeros_like(Plm_est), Nl)
    # alm_l, alm_r = pta.gc2lr_alm(alm_g, alm_c)

    # h_left = pta.inverse_sfft2(alm_l, nside, spin=2)
    # h_right = pta.inverse_sfft2(alm_r, nside, spin=-2)
    # h_plus, h_cross = pta.circ2lin(h_left, h_right) 

    # pta.plot_map(h_plus, name='h-plus', path=path, show=False)
    # pta.plot_map(h_cross, name='h-cross', path=path, show=False)
    # pta.plot_intensities(h_plus, h_cross, path=path, show=False)


    # P_est = pta.inverse_sfft2(Plm_est, nside, spin=0)
    # pta.plot_map(P_est, name='P_est', path=path, show=show)
 
    # # --- Сравнение (опционально) ---
    # diff = np.abs(P_true - P_est)
    # pta.plot_map(diff, name='diff', path=path, show=show)

    # print(f"Максимальная ошибка восстановления: {np.max(diff):.3e}")
    # print(f"Средняя ошибка: {np.mean(diff):.3e}")

    # # # --- Визуализация каталога ---
    # ptacat.plot_catalog(psrcat, path, show=show)

if __name__ == '__main__':
    typer.run(main)