import ptatoolbox as pta
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

dm = pta.DataManager()
path = dm.create_experiment("test-mode")

config = {
    'coords': 'sphere',
    'params': {
        'seed_coord': 42,
        # 'ra_0': 90.0, 
        # 'dec_0': 0.0, 
        # 'alpha': 1.0,
    }
}

# simple_cat = pta.make_catalog(100, {}, mode='simple')
# simple_cat = simple_cat.add(cap_cat.pulsars())
N_psr = 100
simple_cat = pta.make_catalog(N_psr, config, mode='simple')

print(simple_cat)
pta.plot_catalog(simple_cat, path)
pulsars = simple_cat.pulsars()

nside = 100
npix = hp.nside2npix(nside)
vecs = hp.pix2vec(nside, np.arange(npix))
dOmega = hp.nside2pixarea(nside)
Omega = np.column_stack(vecs).T

alpha = - 1
F_plus = alpha * pta.antenna_pattern(pulsars[0], Omega, 1)
F_cross = alpha * pta.antenna_pattern(pulsars[0], Omega, -1)
alm_E, alm_B = hp.map2alm_spin(
    [F_plus, F_cross],
    spin=2
)
alm_e = pta.get_alm_emode(pulsars[0], nside=nside)
direction = - pta.to_vec(pulsars[0])
azimuth, response = pta.average_over_phi(direction, alm_e, nside=nside)
pta.plot_directional_pattern(azimuth, response, path, name='-base')

# pta.plot_map(alm_E, path, name='antenna', catalog=simple_cat, nside=nside)
# pta.plot_map(alm_e, path, name='e-mode', catalog=simple_cat, nside=nside)

# pta.plot_Cl(alm_e, path)
# pta.plot_Cl(alm_E, path)
L = 1000  # максимальное значение l
l_values = np.arange(L + 1)  
counts = 2 * l_values + 1
l_arr = np.repeat(l_values, counts) 
mask = np.logical_and(l_arr != 0, l_arr != 1)
l_arr = l_arr[mask]
Nl = np.sqrt(2/((l_arr+2)*(l_arr+1)*l_arr*(l_arr-1)))

N_modes = 100
modes, vals = pta.eigenmodes(pulsars)
counts = [k+1 for k in range(len(vals))]

num, sigma = counts[:N_modes:], vals[:N_modes:]
sigma_desk = Nl[:len(sigma):] * np.sqrt(N_psr) * (2*np.pi)**0.5
print(sigma[:10:])
print(sigma_desk[:10:])

plt.scatter(num, sigma_desk, color='red', s=20, label='$N_{psr' + '}\\rightarrow \\infty$')
plt.scatter(num, sigma, color='blue', s=20, label='$N_{psr' + '}=' + f'{N_psr}$')
plt.loglog()
plt.legend()
plt.xlabel("mode, m")
plt.ylabel("eigenvalue, $\\sigma_m$")
plt.title("Eigenvalues of PTA modes")
plt.savefig(path / f"EV-{N_psr}.png", dpi=1000)
plt.close()

plus_modes, cross_modes = modes[0], modes[1]
npix = hp.nside2npix(nside)
vecs = hp.pix2vec(nside, np.arange(npix))
dOmega = hp.nside2pixarea(nside)
Omega = np.column_stack(vecs).T
direction = np.array([0, 0, 1])

N = 1
for k in range(N):
    number = (k+1)
    print(number)
    h_plus = plus_modes[k](Omega)
    h_cross = cross_modes[k](Omega)
    sigma = vals[k]
    alm_E, alm_B = hp.map2alm_spin([h_plus, h_cross], spin=2)
    pta.plot_map((-1)*alm_E, path, name=f'e-m{k+1}', catalog=simple_cat)
    # pta.plot_map(alm_e, path, name=f'e-m{k+1}-true')
    # pta.plot_power_mode(plus_modes[k], cross_modes[k], path, catalog=simple_cat)
    # azimuth, response = pta.average_over_phi(direction, alm_E, nside=nside)
    # pta.plot_directional_pattern(azimuth, response, path, name=f'-m{k+1}')
