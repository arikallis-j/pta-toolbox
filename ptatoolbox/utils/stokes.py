import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def map_decomposition(h):
    amplitude = np.abs(h)
    phase = np.mod(np.angle(h), 2 * np.pi)
    return amplitude, phase

def Stokes_parameters(h_plus, h_cross):
    I = np.abs(h_plus)**2 + np.abs(h_cross)**2
    Q = np.abs(h_plus)**2 - np.abs(h_cross)**2
    U = 2 * np.real(h_plus * np.conj(h_cross))
    V = -2 * np.imag(h_plus * np.conj(h_cross))
    return I, Q, U, V

def polarizarion_intensities(h_plus, h_cross):
    I, Q, U, V = Stokes_parameters(h_plus, h_cross)
    I_l = np.sqrt(Q**2 + U**2)
    I_c = np.sqrt(V**2)
    return I_l, I_c

def polarizarion_angles(h_plus, h_cross):
    I, Q, U, V = Stokes_parameters(h_plus, h_cross)
    I_l, I_c = polarizarion_intensities(h_plus, h_cross)
    psi = 0.5 * np.arctan2(U, Q)
    chi = 0.5 * np.arctan2(V, I_l)
    return psi, chi

def polarizarion_degree(h_plus, h_cross):
    I, Q, U, V = Stokes_parameters(h_plus, h_cross)
    P_deg = np.sqrt(np.mean(Q)**2 + np.mean(U)**2 + np.mean(V)**2)/np.mean(I)
    return P_deg

def plot_map(h, name='map', show=True, path=None, grid=True):
    A, Phi = map_decomposition(h)
    plt.figure(figsize=(6, 8))
    hp.mollview(A, sub=211, cmap='viridis', min=0, title='Amplitude')
    if grid:
        hp.graticule()
    hp.mollview(Phi, sub=212, cmap='hsv', min=0, max=2*np.pi, title='Phase')
    if grid:
        hp.graticule()
    plt.tight_layout()
    if path is not None:
        plt.savefig(path / f"{name}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_maps(h, name='maps', show=True, path=None, grid=True):
    R, I = np.real(h), np.imag(h)
    plt.figure(figsize=(6, 8))
    hp.mollview(R, sub=211, cmap='viridis', title='Real')
    if grid:
        hp.graticule()
    hp.mollview(I, sub=212, cmap='viridis', title='Imag')
    if grid:
        hp.graticule()
    plt.tight_layout()
    if path is not None:
        plt.savefig(path / f"{name}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_intensities(h_plus, h_cross, name='stokes-I', show=True, path=None, grid=True, log=False):
    I, Q, U, V = Stokes_parameters(h_plus, h_cross)
    plt.figure(figsize=(6, 8))
    if log:
        norm = LogNorm()
    else:
        norm = None
    hp.mollview(I, sub=211, title='Total Intensity', cmap='inferno', min=0.0, norm=norm)
    if grid:
        hp.graticule()
    hp.mollview(V, sub=212, title='Total Сhirality', cmap='coolwarm', min=-np.max(np.abs(V)), max=np.max(np.abs(V)))
    if grid:
        hp.graticule()
    plt.tight_layout()
    if path is not None:
        plt.savefig(path / f"{name}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_polarization_degrees(h_plus, h_cross, name='stokes-P', show=True, path=None):
    I, Q, U, V = Stokes_parameters(h_plus, h_cross)
    I_l, I_c = polarizarion_intensities(h_plus, h_cross)
    plt.figure(figsize=(6, 8))
    hp.mollview(I_l/I, sub=211, title='Linear Polarization Degree', cmap='plasma', min=0.0, max=1.0)
    hp.graticule()
    hp.mollview(I_c/I, sub=212, title='Circular Polarization Degree', cmap='plasma', min=0.0, max=1.0)
    hp.graticule()
    plt.tight_layout()
    if path is not None:
        plt.savefig(path / f"{name}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_Stokes_parameters(h_plus, h_cross, name='stokes', show=True, path=None):
    I, Q, U, V = Stokes_parameters(h_plus, h_cross)
    psi, chi = polarizarion_angles(h_plus, h_cross)
    plt.figure(figsize=(6, 12))
    hp.mollview(I, sub=311, title='Total Intensity', cmap='inferno', min=0.0)
    hp.graticule()
    hp.mollview(np.mod(2*psi, 2*np.pi), sub=312, title='Linear Polarization Angle $2\psi$', cmap='hsv', min=0, max=2*np.pi)
    hp.graticule()
    hp.mollview(2*chi, sub=313, title='Ellipticity Angle $2\chi$', cmap='RdBu', min=-np.pi/2, max=+np.pi/2)
    hp.graticule()
    plt.tight_layout()
    if path is not None:
        plt.savefig(path / f"{name}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def angular_spectrum(f_lm, nside):
    L = 2 * nside
    cl = np.zeros(L, dtype=np.float64)
    idx = 0
    for l in range(L):
        n_m = 2 * l + 1
        alm_l = f_lm[idx : idx + n_m]
        cl[l] = np.mean(np.abs(alm_l) ** 2)
        idx += n_m
    return cl 

def plot_angular_spectrum(cl, L, log=False, lcut=None):
    l_range = np.arange(L)
    cl_norm = cl

    if lcut is not None:
        l_range = l_range[:lcut:]
        cl_norm = cl_norm[:lcut:]
    
    plt.plot(l_range, cl_norm)
    plt.xlabel('l')
    plt.ylabel('$C_l$')
    if log:
        plt.loglog()
    plt.grid(True)
    plt.show()

def plot_Cl(cl, path = None, name='cl', l_start=2, l_max=1000, show=False):
    ell = np.arange(len(cl))
    l_mean = np.sum(cl * ell) / np.sum(cl)
    f_peak = np.max(cl) / np.sum(cl)
    print(f"l_mean = {l_mean:.2f}")
    print(f"f_mean = {f_peak:.2f}")

    el = ell[l_start:l_max:]
    Nl = np.sqrt(2/((el+2)*(el+1)*el*(el-1)))
    theory = Nl**2 / 10
    # print(cl[l_start:l_max:]/theory)
    plt.plot(ell[l_start:l_max:], cl[l_start:l_max:], label='E-mode (grad)', color='blue')
    plt.plot(ell[l_start:l_max:], theory,  label='theory',  linestyle='--', color='black')
    
    plt.xlabel('Multipole l')
    plt.ylabel('$C_l$')
    plt.loglog()
    plt.legend()
    
    if path is not None:
        plt.savefig(path / f"{name}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_statistics(h_plus, h_cross, path = None, show = True):
    plot_map(h_plus, name='h-plus', path=path, show=show)
    plot_map(h_cross, name='h-cross', path=path, show=show)
    plot_intensities(h_plus, h_cross, path=path, show=show)
    plot_polarization_degrees(h_plus, h_cross, path=path, show=show)
    plot_Stokes_parameters(h_plus, h_cross, path=path, show=show)