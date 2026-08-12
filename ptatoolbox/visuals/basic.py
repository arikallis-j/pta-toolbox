import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from scipy.special import factorial

def get_cartesian(ra, dec, px):
    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)
    px = np.asarray(px, dtype=float)
    ra = np.deg2rad(ra)
    dec = np.deg2rad(dec)
    D = np.where((px > 0) & np.isfinite(px), 1 / px, np.nan)
    
    x = D * np.cos(dec) * np.cos(ra)
    y = D * np.cos(dec) * np.sin(ra)
    z = D * np.sin(dec)
    return x, y, z

def plot_catalog(catalog, path, show=False, save=True):
    ra, dec = catalog.data['RAJD'], catalog.data['DECJD']
    ra = np.where(ra>=180, ra-360, ra)
    ra *= np.pi/180
    dec *= np.pi/180
    plt.figure(figsize=(8,4))
    plt.subplot(projection="mollweide")
    plt.title(f"{catalog.name} pulsar catalog")
    plt.grid(True)
    plt.plot(ra, dec, 'o', markersize=2)
    if save:
        plt.savefig(path / f"{catalog.name}_cat.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_pulsars(catalog, path, show=False, save=True):
    ra, dec, px = catalog.data['RAJD'], catalog.data['DECJD'], catalog.data['PX']
    x, y, z = get_cartesian(ra, dec, px)
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x, y, z, c='blue', s=20, label='PTA')
    ax.scatter(0, 0, 0, c='red', s=100, label='Observer')

    ax.set_box_aspect([1, 1, 1])

    ax.set_xlabel("x, kpc")
    ax.set_ylabel("y, kpc")
    ax.set_zlabel("z, kpc")
    ax.set_title("Distribution of Pulsar Array")
    ax.legend()
    
    if save:
        plt.savefig(path / f"{catalog.name}_pta.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_ppdot(catalog, path, show=False, save=True):
    f0, f1 = catalog.data['F0'], catalog.data['F1']
    P0 = 1/f0
    P1 = - f1/f0**2
    mask_p0 = np.logical_and(P0<1e2, P0>1e-2)
    mask_p1 = np.logical_and(P1<1e-11, P1>1e-18)
    mask = mask_p0 & mask_p1
    P0_psr = np.exp(np.mean(np.log(P0[mask])))
    P1_psr = np.exp(np.mean(np.log(P1[mask])))
    mask_p0 = np.logical_and(P0<1e-2, P0>1e-3)
    mask_p1 = np.logical_and(P1<1e-18, P1>1e-22)
    mask = mask_p0 & mask_p1
    P0_msc = np.exp(np.mean(np.log(P0[mask])))
    P1_msc = np.exp(np.mean(np.log(P1[mask])))

    print(f"P0: {P0_psr:.2e} | F0: {1/P0_psr:.2e}")
    print(f"P1: {P1_psr:.2e} | F1: {-P1_psr/P0_psr**2:.2e}")
    print(f"P0_msc: {P0_msc:.2e} | F0_msc: {1/P0_msc:.2e}")
    print(f"P1_msc: {P1_msc:.2e}| F1_msc: {-P1_msc/P0_msc**2:.2e}")
    a, b = 2.0, -16.0
    log_x_death = np.linspace(-2.3, 2, 1000)
    log_y_death = a * log_x_death + b
    x_death, y_death = 10.0**(log_x_death), 10.0**(log_y_death)
    plt.figure(figsize=(4, 4))
    plt.subplot()
    plt.title(f"{catalog.name} pulsar catalog")
    plt.grid(True)
    plt.xlabel("P")
    plt.ylabel("P_dot")
    plt.scatter(P0, P1, s=10, alpha=0.1)
    plt.scatter(1.0, 1e-15, color='yellow', s=10)
    plt.scatter(P0_psr, P1_psr, color='red', s=10)
    plt.scatter(P0_msc, P1_msc, color='red', s=10)
    plt.plot(x_death, y_death, color='black')
    plt.xscale('log')
    plt.yscale('log')
    

    if save:
        plt.savefig(path / f"{catalog.name}_ppdot.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_map(alm, path, name='map', catalog = None, nside=50, show=False, save=True, cmap='viridis'):
    Map = hp.alm2map(
        alm,
        nside=nside
    )

    dOmega = hp.nside2pixarea(nside)
    power = np.abs(Map)**2
    norm = np.sum(power * dOmega)
    # print(f"norm = {norm:.2f}")
    # print(f"sigma_true = {np.sqrt(norm)}")
    
    hp.mollview(Map, title=name, cmap=cmap, unit="$E(\Omega)$")
    hp.graticule()

    if catalog is not None:
        ra, dec = catalog.data['RAJD'], catalog.data['DECJD']
        hp.projplot(ra, dec, '*', color='white', markersize=5, lonlat=True)
    if save:
        plt.savefig(path / f"{name}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_Cl(alm, path, l_start=2, l_max=1000, show=False, save=True):
    cl = hp.alm2cl(alm)
    ell = np.arange(len(cl))
    l_mean = np.sum(cl * ell) / np.sum(cl)
    f_peak = np.max(cl) / np.sum(cl)
    print(f"l_mean = {l_mean:.2f}")
    print(f"f_mean = {f_peak:.2f}")

    el = ell[l_start:l_max:]
    Nl = np.sqrt(2/((el+2)*(el+1)*el*(el-1)))
    theory = 4 * np.pi * Nl**2/2
    print(cl[l_start:l_max:]/theory)
    plt.plot(ell[l_start:l_max:], cl[l_start:l_max:], label='E-mode (grad)', color='blue')
    plt.plot(ell[l_start:l_max:], theory,  label='theory',  linestyle='--', color='black')
    
    plt.xlabel('Multipole l')
    plt.ylabel('$C_l$')
    plt.loglog()
    plt.legend()
    
    if save:
        plt.savefig(path / f"Cl.png", dpi=1000)
    if show:
        plt.show()
    plt.close()


def plot_directional_pattern(azimuth, response, path, name="", show=False, save=True):
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8,6))
    ax.plot(azimuth, response/np.max(response), 'b-', linewidth=2)
    ax.fill(azimuth, response/np.max(response), alpha=0.3, color='blue')
    ax.set_theta_zero_location('N')  
    ax.set_theta_direction(-1)
    plt.title('Directional pattern')
    if save:
        plt.savefig(path / f"dp{name}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()


def plot_E_mode(p_mode, c_mode, path, catalog = None, nside=50, show=False, save=True, num=0, cmap='viridis'):
    npix = hp.nside2npix(nside)
    vecs = hp.pix2vec(nside, np.arange(npix))
    dOmega = hp.nside2pixarea(nside)
    Omega_hats = np.column_stack(vecs).T
    h_plus = p_mode(Omega_hats)
    h_cross = c_mode(Omega_hats)

    alm_E, alm_B = hp.map2alm_spin(
        [h_plus, h_cross],
        spin=2
    )

    E_map = hp.alm2map(
        alm_E,
        nside=nside
    )

    print(f"norm = {np.sum(np.abs(E_map)**2 * dOmega):.2f}")
    pol = '+' if num>0 else 'x'

    hp.mollview(E_map, title="$E_{" + f"{np.abs(num)}" + "}$ mode", cmap=cmap, unit="$E(\Omega)$")
    hp.graticule()

    if catalog is not None:
        ra, dec = catalog.data['RAJD'], catalog.data['DECJD']
        hp.projplot(ra, dec, '*', color='white', markersize=5, lonlat=True)
    if save:
        plt.savefig(path / f"e-m{np.abs(num)}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_E_power_mode(p_mode, c_mode, path, catalog = None, nside=50, show=False, save=True, num=0, cmap='viridis'):
    npix = hp.nside2npix(nside)
    vecs = hp.pix2vec(nside, np.arange(npix))
    dOmega = hp.nside2pixarea(nside)
    Omega_hats = np.column_stack(vecs).T
    h_plus = p_mode(Omega_hats)
    h_cross = c_mode(Omega_hats)

    alm_E, alm_B = hp.map2alm_spin(
        [h_plus, h_cross],
        spin=2
    )

    E_map = hp.alm2map(
        alm_E,
        nside=nside
    )

    power = np.abs(E_map)**2
    norm = np.sum(power * dOmega)
    
    print(f"norm = {norm:.2f}")
    pol = '+' if num>0 else 'x'

    hp.mollview(power, title="$S_{" + f"{np.abs(num)}" + "}$ mode", cmap=cmap, unit="$|P|$") # ,max=1.0, min=0.0)
    hp.graticule()

    if catalog is not None:
        ra, dec = catalog.data['RAJD'], catalog.data['DECJD']
        hp.projplot(ra, dec, '*', color='white', markersize=5, lonlat=True)
    if save:
        plt.savefig(path / f"e2-m{np.abs(num)}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()


def plot_power_mode(p_mode, c_mode, path, catalog = None, nside=50, show=False, save=True, num=0, cmap='viridis'):
    npix = hp.nside2npix(nside)
    vecs = hp.pix2vec(nside, np.arange(npix))
    dOmega = hp.nside2pixarea(nside)
    Omega_hats = np.column_stack(vecs).T
    h_plus = p_mode(Omega_hats)
    h_cross = c_mode(Omega_hats)

    alm_E, alm_B = hp.map2alm_spin(
        [h_plus, h_cross],
        spin=2
    )
    alm_B = np.zeros_like(alm_E)

    h_plus, h_cross = hp.alm2map_spin(
        [alm_E, alm_B],
        nside=nside,
        spin=2,
        lmax=3*nside - 1
    )

    power = np.abs(h_plus)**2 + np.abs(h_cross)**2
    
    print(f"norm = {np.sum(power * dOmega):.2f}")
    pol = '+' if num>0 else 'x'

    hp.mollview(power, title="$S_{" + f"{np.abs(num)}" + "}$ mode", cmap=cmap, unit="$|P|$") # ,max=1.0, min=0.0)
    hp.graticule()

    if catalog is not None:
        ra, dec = catalog.data['RAJD'], catalog.data['DECJD']
        hp.projplot(ra, dec, '*', color='white', markersize=5, lonlat=True)
    if save:
        plt.savefig(path / f"h2-m{np.abs(num)}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()


def plot_Cl_mode(p_mode, c_mode, path, nside=50, num=0, l_start=2, l_max=1000, show=False, save=True):
    npix = hp.nside2npix(nside)
    vecs = hp.pix2vec(nside, np.arange(npix))
    dOmega = hp.nside2pixarea(nside)
    Omega_hats = np.column_stack(vecs).T
    p_values = p_mode(Omega_hats)
    c_values = c_mode(Omega_hats)
    print(len(c_values))
    almE, almB = hp.map2alm_spin([p_values, c_values], spin=2)
    print(len(almE))
    cl_EE = hp.alm2cl(almE) # grad-like mode
    cl_BB = hp.alm2cl(almB) # curl-like mode
    cl = cl_EE + cl_BB
    print(len(cl))
    ell = np.arange(len(cl_EE))
    l_mean = np.sum(cl_EE * ell) / np.sum(cl_EE)
    f_peak = np.max(cl_EE) / np.sum(cl_EE)
    print(f"l_mean = {l_mean:.2f}")
    print(f"f_mean = {f_peak:.2f}")

    el = ell[l_start:l_max:]
    Nl = np.sqrt(2/((el+2)*(el+1)*el*(el-1)))
    theory = Nl**2 * 3/2

    plt.plot(ell[l_start:l_max:], cl_EE[l_start:l_max:], label='E-mode (grad)', color='blue')
    plt.plot(ell[l_start:l_max:], cl_BB[l_start:l_max:], label='B-mode (curl)', color='red')
    plt.plot(ell[l_start:l_max:], theory,  label='theory',  linestyle='--', color='black')
    # plt.axvline(ell[np.argmax(cl_EE)], color='black', linestyle='--', label='peak')
    
    plt.xlabel('Multipole l')
    plt.ylabel('$C_l$')
    plt.loglog()
    plt.legend()
    
    if save:
        plt.savefig(path / f"Cl_m{np.abs(num)}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_power_Cl_mode(p_mode, c_mode, path, nside=50, num=0, l_start=0, l_max=100, show=False, save=True):
    npix = hp.nside2npix(nside)
    vecs = hp.pix2vec(nside, np.arange(npix))
    dOmega = hp.nside2pixarea(nside)
    Omega_hats = np.column_stack(vecs).T

    h_plus = p_mode(Omega_hats)
    h_cross = c_mode(Omega_hats)

    alm_E, alm_B = hp.map2alm_spin(
        [h_plus, h_cross],
        spin=2
    )
    alm_B = np.zeros_like(alm_E)

    h_plus, h_cross = hp.alm2map_spin(
        [alm_E, alm_B],
        nside=nside,
        spin=2,
        lmax=3*nside - 1
    )

    power = np.abs(h_plus)**2 + np.abs(h_cross)**2
    alm_P = hp.map2alm(power)
    cl_HH = hp.alm2cl(alm_P)
    ell = np.arange(len(cl_HH))

    plt.plot(ell[l_start:l_max:], cl_HH[l_start:l_max:], label='$|h|^2$', color='black')

    plt.xlabel('Multipole l')
    plt.ylabel('$C_l$')
    plt.loglog()
    plt.legend()
    
    if save:
        plt.savefig(path / f"Cl_h2_m{np.abs(num)}.png", dpi=1000)
    if show:
        plt.show()
    plt.close()