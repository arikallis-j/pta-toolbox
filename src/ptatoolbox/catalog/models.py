import numpy as np

from .io import *
from .pulsar import *
from .distr import *
from .names import *

def get_psr_coord(phi, theta, rho=np.nan):
    ra = np.rad2deg(phi)
    dec = np.rad2deg(np.pi/2.0 - theta)
    px = 1/rho
    return ra, dec, px # (deg, deg, mas)

def get_sph_coord(ra, dec, px=np.nan):
    phi = np.deg2rad(ra)
    theta = np.pi/2.0 - np.deg2rad(dec)
    rho = 1/px
    return phi, theta, rho # (rad, rad, kpc)

def simple_test(n_psr):
    Phi, Theta, Rho = isotropic_sphere(n_psr, 42, np.nan)
    Ra, Dec, Px = get_psr_coord(Phi, Theta, Rho)
    Name = [f"S{k}" for k in range(n_psr)]
    return [
        {'name': name, 'ra': ra, 'dec': dec} 
        for name, ra, dec in zip(Name, Ra, Dec)
    ]

def simple_sphere(n_psr, seed_psr=42, radius=np.nan):
    Phi, Theta, Rho = isotropic_sphere(n_psr, seed_psr, radius)
    Ra, Dec, Px = get_psr_coord(Phi, Theta, Rho)
    Name = get_names(Ra, Dec, prefix='S')
    return [
        {'name': name, 'ra': ra, 'dec': dec, 'px': px} 
        for name, ra, dec, px in zip(Name, Ra, Dec, Px)
    ]

def simple_ball(n_psr, seed_psr=42, radius=1):
    Phi, Theta, Rho = isotropic_ball(n_psr, seed_psr, radius)
    Ra, Dec, Px = get_psr_coord(Phi, Theta, Rho)
    Name = get_names(Ra, Dec, prefix='S')
    return [
        {'name': name, 'ra': ra, 'dec': dec, 'px': px} 
        for name, ra, dec, px in zip(Name, Ra, Dec, Px)
    ]

def simple_cone(n_psr, seed_psr=42, alpha=10.0, ra_0=0.0, dec_0=0.0, radius=1):
    phi_0, theta_0, _ = get_sph_coord(ra_0, dec_0)
    Phi, Theta, Rho = isotropic_cone(n_psr, seed_psr, radius, phi_0, theta_0, np.deg2rad(alpha))
    Ra, Dec, Px = get_psr_coord(Phi, Theta, Rho)
    Name = get_names(Ra, Dec, prefix='S')
    return [
        {'name': name, 'ra': ra, 'dec': dec, 'px': px} 
        for name, ra, dec, px in zip(Name, Ra, Dec, Px)
    ]

def simple_cap(n_psr, seed_psr=42, alpha=10.0, ra_0=0.0, dec_0=0.0, radius=np.nan):
    phi_0, theta_0, _ = get_sph_coord(ra_0, dec_0)
    Phi, Theta, Rho = isotropic_cap(n_psr, seed_psr, radius, phi_0, theta_0, np.deg2rad(alpha))
    Ra, Dec, Px = get_psr_coord(Phi, Theta, Rho)
    Name = get_names(Ra, Dec, prefix='S')
    return [
        {'name': name, 'ra': ra, 'dec': dec, 'px': px} 
        for name, ra, dec, px in zip(Name, Ra, Dec, Px)
    ]

def simple_ring(n_psr, seed_psr=42, alpha=10.0, ra_0=0.0, dec_0=0.0, radius=np.nan):
    phi_0, theta_0, _ = get_sph_coord(ra_0, dec_0)
    Phi, Theta, Rho = isotropic_ring(n_psr, seed_psr, radius, phi_0, theta_0, np.deg2rad(alpha))
    Ra, Dec, Px = get_psr_coord(Phi, Theta, Rho)
    Name = get_names(Ra, Dec, prefix='S')
    return [
        {'name': name, 'ra': ra, 'dec': dec, 'px': px} 
        for name, ra, dec, px in zip(Name, Ra, Dec, Px)
    ]

def make_synthetics(n_psr, method, params):
    pulsars = [] 
    synth_data = methods[method](n_psr, **params)
    for k in range(n_psr):
        psr = Pulsar(**synth_data[k])
        pulsars.append(psr)
    data = make_data(pulsars)
    return data 

methods = {
    'test': simple_test,
    'sphere': simple_sphere,
    'ball': simple_ball,
    'cap': simple_cap,
    'ring': simple_ring,
    'cone': simple_cone,
}