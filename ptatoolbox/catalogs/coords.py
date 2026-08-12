"""Coordinate functions for creating pulsar distributions."""

import numpy as np
from collections import Counter

# Coordinate Systems 

def sph2psr(phi, theta, rho=np.nan):
    ra = np.rad2deg(phi)
    dec = np.rad2deg(np.pi/2.0 - theta)
    px = 1/rho
    return ra, dec, px # (deg, deg, mas)

def psr2sph(ra, dec, px=np.nan):
    phi = np.deg2rad(ra)
    theta = np.pi/2.0 - np.deg2rad(dec)
    rho = 1/px
    return phi, theta, rho # (rad, rad, kpc)

def sph2cart(phi, theta):
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return x, y, z

def sph2vec(phi, theta):
    x, y, z = sph2cart(phi, theta)
    return np.array([x, y, z]).T

def vec2cart(vec):
    x, y, z = vec.T
    return x, y, z

def vec2sph(vec):
    x, y, z = vec2cart(vec)
    r = np.linalg.norm(vec, axis=1)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    phi = np.where(phi < 0, phi + 2*np.pi, phi)
    return phi, theta

def vec2psr(vec):
    phi, theta = vec2sph(vec)
    ra, dec, px = sph2psr(phi, theta)
    return ra, dec

def rotation_matrix_from_z(target):
    z_axis = np.array([0.0, 0.0, 1.0])
    if np.allclose(target, z_axis):
        return np.eye(3)
    axis = np.cross(z_axis, target)
    axis = axis / np.linalg.norm(axis)
    angle = np.arccos(np.dot(z_axis, target))
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)
    return R


# Distributions

def isotropic_sphere(n_array, seed_coord=42, radius=np.nan):
    rng = np.random.default_rng(seed=seed_coord)
    phi = rng.uniform(0, 2*np.pi, n_array)
    cos_theta = rng.uniform(-1, 1, n_array)
    theta = np.arccos(cos_theta)
    rho = radius * np.ones(n_array)
    ra, dec, px = sph2psr(phi, theta, rho)
    return ra, dec, px

def isotropic_ball(n_array, seed_coord=42, radius=np.nan):
    rng = np.random.default_rng(seed=seed_coord)
    phi = rng.uniform(0, 2*np.pi, n_array)
    cos_theta = rng.uniform(-1, 1, n_array)
    theta = np.arccos(cos_theta)
    rho = radius * rng.uniform(0, 1, n_array) ** (1/3)
    ra, dec, px = sph2psr(phi, theta, rho)
    return ra, dec, px

def isotropic_cap(n_array, seed_coord=42, radius=np.nan, ra_0=0.0, dec_0=0.0, alpha=45.0):
    phi_0, theta_0, _ = psr2sph(ra_0, dec_0)
    alpha = np.deg2rad(alpha)
    rng = np.random.default_rng(seed=seed_coord)
    phi_pole = rng.uniform(0, 2*np.pi, n_array)
    u_pole = rng.uniform(0, 1, n_array)
    theta_pole = np.arccos(1 - u_pole * (1 - np.cos(alpha)))
    
    points_pole = sph2vec(phi_pole, theta_pole)
    axis = sph2vec(phi_0, theta_0)
    R = rotation_matrix_from_z(axis)
    points = points_pole @ R.T
    phi, theta = vec2sph(points)
    rho = radius * np.ones(n_array)
    ra, dec, px = sph2psr(phi, theta, rho)
    return ra, dec, px

def isotropic_cone(n_array, seed_coord=42, radius=np.nan, ra_0=0.0, dec_0=0.0, alpha=45.0):
    phi_0, theta_0, _ = psr2sph(ra_0, dec_0)
    alpha = np.deg2rad(alpha)
    rng = np.random.default_rng(seed=seed_coord)
    phi_pole = rng.uniform(0, 2*np.pi, n_array)
    u_pole = rng.uniform(0, 1, n_array)
    theta_pole = np.arccos(1 - u_pole * (1 - np.cos(alpha)))
    
    points_pole = sph2vec(phi_pole, theta_pole)
    axis = sph2vec(phi_0, theta_0)
    R = rotation_matrix_from_z(axis)
    points = points_pole @ R.T
    phi, theta = vec2sph(points)
    rho = radius * rng.uniform(0, 1, n_array) ** (1/3)
    ra, dec, px = sph2psr(phi, theta, rho)
    return ra, dec, px

def isotropic_ring(n_array, seed_coord=42, radius=np.nan, ra_0=0.0, dec_0=0.0, alpha=45.0):
    phi_0, theta_0, _ = psr2sph(ra_0, dec_0)
    alpha = np.deg2rad(alpha)
    rng = np.random.default_rng(seed=seed_coord)
    phi_pole = rng.uniform(0, 2*np.pi, n_array)
    u_pole = np.ones(n_array)
    theta_pole = np.arccos(1 - u_pole * (1 - np.cos(alpha)))
    
    points_pole = sph2vec(phi_pole, theta_pole)
    axis = sph2vec(phi_0, theta_0)
    R = rotation_matrix_from_z(axis)
    points = points_pole @ R.T
    phi, theta = vec2sph(points)
    rho = radius * np.ones(n_array)
    ra, dec, px = sph2psr(phi, theta, rho)
    return ra, dec, px

COORDS = {
    'sphere': {
        'distr': isotropic_sphere,
        'params': ['seed_coord', 'radius'],
    },
    'ball': {
        'distr': isotropic_ball,
        'params': ['seed_coord', 'radius'],
    },
    'cap':{
        'distr': isotropic_cap,
        'params': ['seed_coord', 'radius', 'ra_0', 'dec_0', 'alpha']
    }, 
    'cone': {
        'distr': isotropic_cone,
        'params': ['seed_coord', 'radius', 'ra_0', 'dec_0', 'alpha'],
    },
    'ring': {
        'distr': isotropic_ring,
        'params': ['seed_coord', 'radius', 'ra_0', 'dec_0', 'alpha'],
    },
}