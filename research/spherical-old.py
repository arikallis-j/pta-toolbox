# Импорты
import jax
jax.config.update("jax_enable_x64", True)  # Для научной точности
import jax.numpy as jnp
import numpy as np
import s2fft
import healpy as hp
import matplotlib.pyplot as plt

L = 2

F_lm = np.zeros((L, 2 * L - 1), dtype=np.complex128)
for l in range(L):
    for m in range(-l, l + 1):
        idx = L - 1 + m
        F_lm[l, idx] = l+1j*m
print(F_lm)

f_lm = s2fft.sampling.s2_samples.flm_2d_to_1d(F_lm, L)
print(f_lm)

L1 = 3
g_lm = np.zeros((L1**2), dtype=np.complex128)
g_lm[:L**2] = f_lm
print(g_lm)

G_lm = s2fft.sampling.s2_samples.flm_1d_to_2d(g_lm, L1)
print(G_lm)

L2 = 4
h_lm = np.zeros((L2**2), dtype=np.complex128)
for l in range(L2):
    for m in range(-l, l + 1):
        idx = l*(l+1) + m
        h_lm[idx] = l+1j*m
print(h_lm)

def make_map(f_lm, nside=10, spin=0):
    Lmax = 3*nside - 1
    map_lm = np.zeros((Lmax**2), dtype=np.complex128)
    map_lm[:len(f_lm):] = f_lm
    Map_lm = s2fft.sampling.s2_samples.flm_1d_to_2d(map_lm, Lmax)
    map_alpha = s2fft.inverse_jax(
        Map_lm, Lmax, 
        sampling="healpix", 
        nside=nside, 
        spin=spin, 
        L_lower=abs(spin)
    )
    map_alpha = np.array(map_alpha)
    return map_alpha

def make_background(L):
    ell = np.arange(L)
    C_l = np.ones((L), dtype=np.complex128)
    N_l = np.zeros((L), dtype=np.complex128)
    N_l[ell>1] = np.sqrt(2/((ell[ell>1]-1)*ell[ell>1]*(ell[ell>1]+1)*(ell[ell>1]+2)))
    a_lm = np.zeros((L**2), dtype=np.complex128)
    psi_lm = np.zeros((L**2), dtype=np.complex128)
    for l in range(L):
        for m in range(-l, l + 1):
            idx = l*(l+1) + m
            std = np.sqrt(C_l[l])
            real = np.random.normal(0.0, std)
            imag = np.random.normal(0.0, std)
            a_lm[idx] = real + 1j * imag
            psi_lm[idx] = a_lm[idx] * N_l[l] / 2

    return a_lm, psi_lm

L2 = 2

g_lm, psi_g_lm = make_background(L2)
c_lm, psi_c_lm = make_background(L2)

a_L_lm = g_lm #+ 1j*c_lm
a_R_lm = g_lm #- 1j*c_lm

print(g_lm)
print(c_lm)

nside = 50
map_a = make_map(g_lm, nside=nside)
map_psi = make_map(psi_g_lm, nside=nside)
map_L = make_map(a_L_lm, nside=nside, spin=+2)
map_R = make_map(a_R_lm, nside=nside, spin=-2) 
map_P = np.abs(map_R)**2 + np.abs(map_L)**2

MAP = map_P

print(f"Карта на сетке HEALPix имеет {len(MAP)} пикселей")
# 3. Визуализация амплитуды и фазы
# Амплитуда
amplitude = np.abs(MAP)
# Фаза (в радианах, от -π до π)
phase = np.angle(MAP)

# 1. Приводим фазу в диапазон [0, 2π)
phase_shifted = np.mod(phase, 2 * np.pi)
# 2. Находим максимум амплитуды для задания верхней границы
amp_max = np.max(amplitude)

# 3. Визуализация
plt.figure(figsize=(12, 6))

# --- Амплитуда ---
plt.subplot(1, 2, 1)
hp.mollview(
    amplitude,
    title="Амплитуда (|ψ(α)|)",
    cmap="viridis",
    min=0,               # нижняя граница
    max=amp_max,         # верхняя граница = максимум
    hold=True
)
plt.gca().set_title("Амплитуда", fontsize=14)

# --- Фаза ---
plt.subplot(1, 2, 2)
hp.mollview(
    phase_shifted,
    title="Фаза (arg ψ(α))",
    cmap="hsv",          # циклическая карта для фазы
    min=0,               # от 0
    max=2 * np.pi,       # до 2π
    hold=True
)
plt.gca().set_title("Фаза", fontsize=14)

plt.tight_layout()
plt.show()