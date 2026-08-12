import numpy as np
import healpy as hp

from scipy.special import sph_harm_y
from tqdm import tqdm

def to_unit(ra, dec):
    return np.array([
        np.cos(dec) * np.cos(ra),
        np.cos(dec) * np.sin(ra),
        np.sin(dec)
    ])

def to_vec(psr):
    u = to_unit(np.deg2rad(psr.ra), np.deg2rad(psr.dec))
    return u
    
def angular_distance(u1, u2):
    cos_theta = np.clip(np.dot(u1, u2), -1.0, 1.0)
    return np.arccos(cos_theta)

def isotropic_hd(p1, p2):
    x = np.clip(np.dot(p1, p2), -1.0, 1.0)
    if np.isclose(angular_distance(p1, p2), 0.0):
        return 1.0
    return 1.0 + 3/2 * (1 - x) * (np.log((1 - x)/2) - 1/6)

def calc_ORF(pulsars, hellings_downs, params):
    Npsr = len(pulsars)
    Gamma = np.zeros((Npsr, Npsr), dtype=float)
    for i in tqdm(range(Npsr)):
        for j in range(Npsr):
            p1, p2 = to_vec(pulsars[i]), to_vec(pulsars[j])
            Gamma[i, j] = hellings_downs(p1, p2, **params)
    return Gamma

import numpy as np

def get_polarization_basis(Omega):
    Omega = np.asarray(Omega)

    single = False
    if Omega.ndim == 1:
        Omega = Omega[:, None]
        single = True

    x, y, z = Omega[0], Omega[1], Omega[2]

    phi = np.arctan2(y, x)
    theta = np.arccos(np.clip(z, -1.0, 1.0))

    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)

    e_theta = np.vstack((
        cos_theta * np.cos(phi),
        cos_theta * np.sin(phi),
        -sin_theta
    ))

    e_phi = np.vstack((
        -np.sin(phi),
        np.cos(phi),
        np.zeros_like(phi)
    ))

    if single:
        return e_theta[:, 0], e_phi[:, 0]

    return e_theta, e_phi
    
def antenna_pattern(pulsar, Omega, A):
    p = to_vec(pulsar)
    u, v = get_polarization_basis(Omega)
    p_dot_u = np.dot(p, u)
    p_dot_v = np.dot(p, v)
    Omega_dot_p = np.dot(p, Omega)

    denom = 1.0 + Omega_dot_p
    F_plus = 1/2 * (p_dot_u**2 - p_dot_v**2) / denom
    F_cross = 1/2 * (2 * p_dot_u * p_dot_v) / denom
    
    if A == 1:
        return F_plus
    elif A == -1:
        return F_cross
    else:
        return None

def antenna_correlation(psr_1, psr_2, nside=10):
    npix = hp.nside2npix(nside)
    vecs = hp.pix2vec(nside, np.arange(npix))
    dOmega = hp.nside2pixarea(nside)
    Omega = np.column_stack(vecs).T
    M = 0.0
    for a in [-1, 1]:
        F_1 = antenna_pattern(psr_1, Omega, a)
        F_2 = antenna_pattern(psr_2, Omega, a)
        M += np.sum(F_1 * F_2 * dOmega)
    return M

def calc_antenna_Fisher(pulsars, nside=10):
    Npsr = len(pulsars)
    Gamma = np.zeros((Npsr, Npsr), dtype=float)
    for i in tqdm(range(Npsr)):
        for j in range(Npsr):
            p1, p2 = pulsars[i], pulsars[j]
            Gamma[i, j] = antenna_correlation(p1, p2, nside)
    return Gamma

def eigenmode(pulsars, Omega, A, s_val, u_vec):
    V_map, norm = 0, s_val
    for i in range(len(pulsars)):
        psr, amp = pulsars[i], u_vec[i]
        V_map += amp * antenna_pattern(psr, Omega, A)
    return V_map / norm

def eigenmodes(pulsars):
    M = (4*np.pi)/3 * calc_ORF(pulsars, isotropic_hd, {})
    # M = calc_antenna_Fisher(pulsars)
    eigenval, eigenvec = np.linalg.eigh(M)

    # сортировка по убыванию
    idx = np.argsort(eigenval)[::-1]
    eigenval = eigenval[idx]
    eigenvec = eigenvec[:, idx]
    sigmaval = np.sqrt(eigenval)

    pols = [1, -1]
    modes = []
    for a in range(len(pols)):
        a_modes = []
        for k in range(len(eigenval)):
            A_pol, s_val, u_vec = pols[a], sigmaval[k], eigenvec[:,k]
            def v_ak(Omega, A_pol=A_pol, s_val=s_val, u_vec=u_vec):
                return eigenmode(pulsars, Omega, A_pol, s_val, u_vec)
            a_modes.append(v_ak)
        modes.append(a_modes)
    return modes, sigmaval
        

def get_EB_from_polarization(h_plus, h_cross, nside):
    # spin-2 spherical transform
    alm_E, alm_B = hp.map2alm_spin(
        [h_plus, h_cross],
        spin=2
    )

    # восстановление карт E и B
    E_map = hp.alm2map(
        alm_E,
        nside=nside
    )

    B_map = hp.alm2map(
        alm_B,
        nside=nside
    )

    return E_map, B_map, alm_E, alm_B

def E_map_to_polarization(E_map, nside=None, lmax=None):
    if nside is None:
        nside = hp.get_nside(E_map)

    if lmax is None:
        lmax = 3*nside - 1

    alm_E = hp.map2alm(
        E_map,
        lmax=lmax
    )

    # 2. Чистая E-мода:
    alm_B = np.zeros_like(alm_E)

    # 3. Обратное spin-2 преобразование
    h_plus, h_cross = hp.alm2map_spin(
        [alm_E, alm_B],
        nside=nside,
        spin=2,
        lmax=lmax
    )

    return h_plus, h_cross, alm_E

def alm_delta(l, m, psr):
    vec = to_vec(psr)
    theta, phi = hp.vec2ang(vec)
    Nl = 1.0
    alm = 4 * np.pi * (-1)**l * Nl/np.sqrt(2) * np.conj(sph_harm_y(l, m, theta, phi))
    return alm

def alm_gauss(l, m, psr):
    vec = to_vec(psr)
    theta, phi = hp.vec2ang(vec)
    Nl = np.sqrt(2/((l+2)*(l+1)*l*(l-1)))

    Norm = 4 * np.pi/np.sqrt(2)
    alm = Norm * (-1)**l * Nl * np.conj(sph_harm_y(l, m, theta, phi))
    return alm


def get_alm_emode(psr, nside=100, lmax=None):
    if lmax is None:
        lmax = 3 * nside - 1

    n_alm = hp.Alm.getsize(lmax)
    alms = np.zeros(n_alm, dtype=np.complex128)
    vals = np.zeros(n_alm, dtype=np.complex128)
    l_arr = np.repeat(np.arange(lmax + 1), np.arange(1, lmax + 2))
    m_arr = np.concatenate([np.arange(l + 1) for l in range(lmax + 1)])
    mask = np.logical_and(l_arr != 0, l_arr != 1)
    idxs = hp.Alm.getidx(lmax, l_arr, m_arr)

    vals[mask] = alm_gauss(l_arr[mask], m_arr[mask], psr)
    vals[~mask] = 0.0
    alms[idxs] = vals
    return alms

def average_over_phi_cos(psr, alm, nside=100, nbins=100):
    npix = hp.nside2npix(nside)
    vecs = hp.pix2vec(nside, np.arange(npix))
    Omega = np.column_stack(vecs).T

    E_map = hp.alm2map(alm, nside=nside)
    p = - to_vec(psr)

    cos_theta = np.dot(p, Omega)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    
    bins = np.linspace(-1.0, 1.0, nbins+1)
    bin_indices = np.digitize(cos_theta, bins) - 1
    
    sum_vals, _ = np.histogram(cos_theta, bins=bins, weights=E_map)
    count, _ = np.histogram(cos_theta, bins=bins)
    mask = count > 0
    E_mean = np.zeros(nbins)
    E_mean[mask] = sum_vals[mask] / count[mask]
    
    cos_centers = (bins[:-1] + bins[1:]) / 2.0
    theta_centers = np.arccos(cos_centers)

    north_pole = np.array([0.0, 0.0, 1.0])
    south_pole = np.array([0.0, 0.0, -1.0])

    # Индексы пикселей, содержащих эти направления
    ipix_north = hp.vec2pix(nside, *north_pole)
    ipix_south = hp.vec2pix(nside, *south_pole)

    # Значения карты в этих пикселях
    E_north = E_map[ipix_north]
    E_south = E_map[ipix_south]

    theta_full = np.concatenate([[0.0], theta_centers[::-1], [np.pi], 2*np.pi - theta_centers])   # исключаем дублирование нуля
    E_full = np.concatenate([[E_north], E_mean[::-1], [E_south], E_mean]) # зеркалим значения
    
    return theta_full, E_full

def average_over_phi(direction, alm, nside=100, nbins=100):
    npix = hp.nside2npix(nside)
    vecs = hp.pix2vec(nside, np.arange(npix))
    Omega = np.column_stack(vecs).T

    E_map = np.abs(hp.alm2map(alm, nside=nside))
    cos_theta = np.dot(direction, Omega)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)

    # Биннинг по θ (равномерные бины от 0 до π)
    bins = np.linspace(0.0, np.pi, nbins + 1)
    bin_indices = np.digitize(theta, bins) - 1   # индексы от 0 до nbins-1
    
    # Суммируем значения и считаем количество пикселей в каждом бине
    sum_vals = np.zeros(nbins)
    count = np.zeros(nbins, dtype=int)
    for i in range(nbins):
        mask = (bin_indices == i)
        if np.any(mask):
            sum_vals[i] = np.sum(E_map[mask])
            count[i] = np.sum(mask)
    
    # Вычисляем среднее (там, где есть пиксели)
    valid = count > 0
    E_mean = np.zeros(nbins)
    E_mean[valid] = sum_vals[valid] / count[valid]
    
    # Центры бинов по θ
    theta_centers = (bins[:-1] + bins[1:]) / 2.0

    north_pole = + direction
    south_pole = - direction

    # Индексы пикселей, содержащих эти направления
    ipix_north = hp.vec2pix(nside, *north_pole)
    ipix_south = hp.vec2pix(nside, *south_pole)

    # Значения карты в этих пикселях
    E_north = E_map[ipix_north]
    E_south = E_map[ipix_south]

    theta_full = np.concatenate([theta_centers[::-1], [0.0], 2*np.pi - theta_centers,])   # исключаем дублирование нуля
    E_full = np.concatenate([E_mean[::-1], [E_north], E_mean]) # зеркалим значения
    
    return theta_full, E_full