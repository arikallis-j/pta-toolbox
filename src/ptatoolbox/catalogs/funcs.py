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

def isotropic_ball(n_array, seed, radius):
    rng = np.random.default_rng(seed=seed)
    phi = rng.uniform(0, 2*np.pi, n_array)
    cos_theta = rng.uniform(-1, 1, n_array)
    theta = np.arccos(cos_theta)
    rho = radius * rng.uniform(0, 1, n_array) ** (1/3)
    return phi, theta, rho

def isotropic_sphere(n_array, seed, radius):
    rng = np.random.default_rng(seed=seed)
    phi = rng.uniform(0, 2*np.pi, n_array)
    cos_theta = rng.uniform(-1, 1, n_array)
    theta = np.arccos(cos_theta)
    rho = radius * np.ones(n_array)
    return phi, theta, rho

def isotropic_cone(n_array, seed, radius, phi_0, theta_0, alpha):
    rng = np.random.default_rng(seed=seed)
    phi_pole = rng.uniform(0, 2*np.pi, n_array)
    u_pole = rng.uniform(0, 1, n_array)
    theta_pole = np.arccos(1 - u_pole * (1 - np.cos(alpha)))
    
    points_pole = sph2vec(phi_pole, theta_pole)
    axis = sph2vec(phi_0, theta_0)
    R = rotation_matrix_from_z(axis)
    points = points_pole @ R.T
    phi, theta = vec2sph(points)
    rho = radius * rng.uniform(0, 1, n_array) ** (1/3)
    return phi, theta, rho

def isotropic_cap(n_array, seed, radius, phi_0, theta_0, alpha):
    rng = np.random.default_rng(seed=seed)
    phi_pole = rng.uniform(0, 2*np.pi, n_array)
    u_pole = rng.uniform(0, 1, n_array)
    theta_pole = np.arccos(1 - u_pole * (1 - np.cos(alpha)))
    
    points_pole = sph2vec(phi_pole, theta_pole)
    axis = sph2vec(phi_0, theta_0)
    R = rotation_matrix_from_z(axis)
    points = points_pole @ R.T
    phi, theta = vec2sph(points)
    rho = radius * np.ones(n_array)
    
    return phi, theta, rho

def isotropic_ring(n_array, seed, radius, phi_0, theta_0, alpha):
    rng = np.random.default_rng(seed=seed)
    phi_pole = rng.uniform(0, 2*np.pi, n_array)
    u_pole = np.ones(n_array)
    theta_pole = np.arccos(1 - u_pole * (1 - np.cos(alpha)))
    
    points_pole = sph2vec(phi_pole, theta_pole)
    axis = sph2vec(phi_0, theta_0)
    R = rotation_matrix_from_z(axis)
    points = points_pole @ R.T
    phi, theta = vec2sph(points)
    rho = radius * np.ones(n_array)
    
    return phi, theta, rho

# Pulsar Names

def get_name(ra, dec, prefix='S'):
    ra_h = np.floor(ra/15.0).astype(int)
    ra_m = np.floor(np.round((ra/15.0 - ra_h) * 60, decimals=1)).astype(int)
    dec_d = np.floor(np.abs(dec)).astype(int)
    dec_m = np.floor(np.round((np.abs(dec) - dec_d) * 60, decimals=1)).astype(int)
    dec_s = "+" if np.sign(dec)>=0.0 else "-"
    name = f"{prefix}{ra_h:02d}{ra_m:02d}{dec_s}{dec_d:02d}{dec_m:02d}"
    return name

def get_names(ra, dec, prefix='S'):
    base_names = [get_name(r, d, prefix) for r, d in zip(ra, dec)]
    freq = Counter(base_names)
    counters = {name: 0 for name in freq}
    final_names = []
    for name in base_names:
        cnt = counters[name]
        if cnt == 0:
            final_names.append(name)
        else:
            if cnt <= 26:
                suffix = chr(ord('A') + cnt - 1)
            else:
                idx = cnt - 27
                first = chr(ord('a') + idx // 26)
                second = chr(ord('a') + idx % 26)
                suffix = first + second
            final_names.append(name + suffix)
        counters[name] += 1
    
    return final_names