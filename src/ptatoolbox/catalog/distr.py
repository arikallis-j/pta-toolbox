import numpy as np

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