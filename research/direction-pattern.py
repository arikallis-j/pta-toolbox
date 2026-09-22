import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.special import spherical_jn, eval_legendre

import ptatoolbox as pta 
dm = pta.DataManager()

def unit_vector_from_angles(theta, phi):
    """
    Вектор единичной длины по сферическим углам.
    theta: полярный угол [0, pi]
    phi: азимут [0, 2pi)
    """
    return np.array([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta)
    ], dtype=float)


def normalize(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("Нельзя нормализовать нулевой вектор.")
    return v / n


def F_l_G(l, y):
    """
    \mathcal{F}_l^G(y)

    Важно:
    - scipy.special.spherical_jn не умеет j_{-1}, поэтому для l=0
      используем аналитическое продолжение j_{-1}(y) = cos(y)/y.
    - Формула корректна при y != 0. Для y=0 возникает особенность из-за 1/y.
    """
    y = complex(y)

    if y == 0:
        raise ValueError("Для этой формулы нужно y != 0.")

    jl = spherical_jn(l, y)

    if l == 0 or l == 1:
        return 0.0
    else:
        j_l_minus_1 = spherical_jn(l - 1, y)

    pref = (1j ** l) * np.exp(-1j * y)

    term1 = (y - 1j * l * (l + 1) / 2.0) * j_l_minus_1
    term2 = (
        (l * (l - 1) / 2.0)
        + 1j * l * (l**2 - 1) / (2.0 * y)
        + 1j * y
    ) * jl

    return pref * (term1 + term2)


def K_on_vectors(y, khat, khat_prime, lmax):
    """
    K(y, khat, khat') для массива направлений khat.

    Parameters
    ----------
    y : float
    khat : array, shape (..., 3)
        Массив единичных векторов на сфере.
    khat_prime : array, shape (3,)
        Константный единичный вектор.
    lmax : int
        Максимальное l в сумме.

    Returns
    -------
    K : complex ndarray, shape (...)
    """
    khat = np.asarray(khat, dtype=float)
    khat_prime = normalize(khat_prime)

    if khat.shape[-1] != 3:
        raise ValueError("khat должен иметь форму (..., 3).")

    # cos(gamma) = k · k'
    mu = np.tensordot(khat, khat_prime, axes=([-1], [0]))

    out = np.zeros(mu.shape, dtype=np.complex128)

    for l in range(lmax + 1):
        Fl = F_l_G(l, y)
        Pl = eval_legendre(l, mu)  # P_l(mu), работает и с массивом mu
        out += (2 * l + 1) / (4.0 * np.pi) * Fl * Pl

    return out


def K_healpix_map(y, khat_prime, nside, lmax, nest=False):
    """
    Считает K на всей HEALPix-карте.
    """
    npix = hp.nside2npix(nside)
    pix = np.arange(npix)

    x, yy, z = hp.pix2vec(nside, pix, nest=nest)
    khat = np.stack([x, yy, z], axis=-1)

    return K_on_vectors(y, khat, khat_prime, lmax)

def K_vs_theta(y, theta, lmax):
    """
    Вычисляет K(y, θ) для угла θ ∈ [0, π] между \hat k и \hat k'.
    μ = cos θ.
    Возвращает комплексный массив той же формы, что и theta.
    """
    theta = np.asarray(theta)
    mu = np.cos(theta)
    out = np.zeros_like(mu, dtype=np.complex128)
    for l in range(lmax + 1):
        Fl = F_l_G(l, y)
        Pl = eval_legendre(l, mu)          # Pl(μ)
        out += (2 * l + 1) / (4.0 * np.pi) * Fl * Pl
        print(f"{l}/{lmax}", np.max(np.abs(Fl)))
    return out


def plot_directional_pattern(path, y, lmax, npoints=500, save_dir="data", with_phase=True):
    """
    Строит диаграмму направленности:
      - левый подграфик: |K|(θ) в полярных координатах,
      - правый подграфик: arg(K)(θ) в полярных координатах (радиус = 1, цвет = фаза).
    Сохраняет рисунок в save_dir/directional_pattern.png.
    """
    theta = np.linspace(0, 2*np.pi, npoints)
    Kvals = K_vs_theta(y, theta, lmax)
    amp = np.abs(Kvals)
    phase = np.angle(Kvals)                # [-π, π]
    phase = np.mod(phase, 2.0 * np.pi)     # [0, 2π)

    fig, ax1 = plt.subplots(
        1, 1, subplot_kw={'projection': 'polar'}, figsize=(7, 6)
    )

    r_phase = amp
    ax1.fill(theta, amp, alpha=0.3, color='blue')
    if with_phase:
        sc = ax1.scatter(theta, r_phase, c=phase, cmap='hsv',
                     vmin=0, vmax=2*np.pi, s=20, alpha=0.9)

        # Цветовая шкала для фазы
        cbar = fig.colorbar(sc, ax=ax1, orientation='vertical', pad=0.1)
        cbar.set_label('phase [rad]')
    else:
        ax1.plot(theta, amp, 'b-', linewidth=1)
    ax1.set_title(r'Directional Pattern  $K(\theta)$', fontsize=14)
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)
    ax1.grid(True)

    plt.tight_layout()
    plt.savefig(f"{path}/DP-{y:.0f}y.png", dpi=300, bbox_inches='tight')
    plt.show()


def plot_directional_filter(path, y, lmax, npoints=500, save_dir="data", with_phase=True):
    """
    Строит диаграмму направленности:
      - левый подграфик: |K|(θ) в полярных координатах,
      - правый подграфик: arg(K)(θ) в полярных координатах (радиус = 1, цвет = фаза).
    Сохраняет рисунок в save_dir/directional_pattern.png.
    """
    theta = np.linspace(-np.pi, np.pi, npoints)
    Kvals = K_vs_theta(y, theta, lmax)
    amp = np.abs(Kvals)
    phase = np.angle(Kvals)                # [-π, π]
    phase = np.mod(phase, 2.0 * np.pi)     # [0, 2π)
    theta_max = np.abs(theta[np.argmax(amp)])
    print(f"angle = {theta_max * 180/np.pi:.4f} deg")
    fig, ax1 = plt.subplots(1, 1, figsize=(7, 6))

    r_phase = amp
    if with_phase:
        sc = ax1.scatter(theta, r_phase, c=phase, cmap='hsv',
                     vmin=0, vmax=2*np.pi, s=20, alpha=0.9)
        # Цветовая шкала для фазы
        cbar = fig.colorbar(sc, ax=ax1, orientation='vertical', pad=0.1)
        cbar.set_label('phase [rad]')
    else:
        ax1.plot(theta, amp, 'b-', linewidth=1)
    ax1.axvline(+theta_max, color='black', linestyle='--')
    ax1.axvline(-theta_max, color='black', linestyle='--')
    ax1.set_title(r'Directional Filter  $K(\theta)$', fontsize=14)
    ax1.grid(True)

    plt.tight_layout()
    plt.savefig(f"{path}/DF-{lmax}l.png", dpi=300, bbox_inches='tight') # {y:.0f}y
    # plt.show()

def plot_K_amplitude(path, y, khat_prime, nside, lmax, nest=False, cmap="viridis"):
    """
    Рисует карту |K|.
    """
    Kmap = K_healpix_map(y, khat_prime, nside, lmax, nest=nest)
    amp = np.abs(Kmap)

    # plt.figure(figsize=(10, 6))
    hp.mollview(
        amp,
        title=rf"$|K|,\ y={y},\ nside={nside},\ lmax={lmax}$",
        cmap=cmap,
        nest=nest,
        cbar=True,
    )
    hp.graticule()
    plt.savefig(f"{path}/Amp.png")

    return Kmap, amp


def plot_K_phase(path, y, khat_prime, nside, lmax, nest=False):
    """
    Рисует карту arg(K) в диапазоне [0, 2pi).
    """
    Kmap = K_healpix_map(y, khat_prime, nside, lmax, nest=nest)

    phase = np.mod(np.angle(Kmap), 2.0 * np.pi)

    # Можно слегка замаскировать точки, где амплитуда почти нулевая,
    # потому что фаза там численно плохо определена.
    amp = np.abs(Kmap)
    phase = np.ma.array(phase, mask=(amp < amp.max() * 1e-12))

    # plt.figure(figsize=(10, 6))
    hp.mollview(
        phase,
        title=rf"$\arg K \in [0,2\pi),\ y={y},\ nside={nside},\ lmax={lmax}$",
        cmap="hsv",
        min=0.0,
        max=2.0 * np.pi,
        nest=nest,
        cbar=True,
    )
    hp.graticule()
    plt.savefig(f"{path}/Phase.png")

    return Kmap, phase


AU = 150_000_000
PC = 206265 * AU
KPC = 1000 * PC
C_LIGHT = 300_000
F_NHZ = 10**(-9)
PI = np.pi

L = 1 # kpc
f = 1 # nHz
omega = 2*PI*f
y = omega * L * KPC*F_NHZ/C_LIGHT
print(y)

theta_p = np.deg2rad(90.0)
phi_p = np.deg2rad(0.0)
khat_prime = unit_vector_from_angles(theta_p, phi_p)
path = dm.create_experiment("DF-648y-ph")

# nside = 1000
for lmax in range(2, 50):#int(y*1.1) 
    plot_directional_filter(path, y, lmax, npoints=10_000, with_phase=True)
# plot_directional_pattern(path, y, lmax, npoints=1_000, with_phase=True)

# Kmap, amp = plot_K_amplitude(path, y, khat_prime, nside, lmax)
# Kmap, phase = plot_K_phase(path, y, khat_prime, nside, lmax)
