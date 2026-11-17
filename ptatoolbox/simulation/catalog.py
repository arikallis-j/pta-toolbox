import numpy as np
import pta_catalog as ptacat

def make_catalog(n_psr, config, mode='simple'):
    catalog = ptacat.make_catalog(n_psr, config, mode=mode)
    return catalog

def make_pulsars(n_psr, config, mode='simple'):
    catalog = ptacat.make_catalog(n_psr, config, mode=mode)
    return catalog.data

def make_coords(n_psr, config, mode='simple'):
    catalog = ptacat.make_catalog(n_psr, config, mode=mode)
    data = catalog.data
    RA, DEC, PX = data['RAJD'].values, data['DECJD'].values, data['PX'].values
    PHI, THETA, RHO = ptacat.psr2sph(RA, DEC, PX) # (rad, rad, kpc)
    return PHI, THETA, RHO 

def get_coords(data):
    RA, DEC, PX = data['RAJD'].values, data['DECJD'].values, data['PX'].values
    PHI, THETA, RHO = ptacat.psr2sph(RA, DEC, PX) # (rad, rad, kpc)
    return PHI, THETA, RHO 

SPHERE = {
    'coords': 'sphere',
    'params': {
        'seed_coord': 42,
        'radius': np.nan,
    }
}

BALL = {
    'coords': 'ball',
    'params': {
        'seed_coord': 42,
        'radius': 1.0,
    }
}

CAP = {
    'coords': 'cap',
    'params': {
        'seed_coord': 42,
        'radius': np.nan,
        'ra_0': 0.0,
        'dec_0': 0.0, 
        'alpha': 30.0,
    }
}

CONE = {
    'coords': 'cone',
    'params': {
        'seed_coord': 42,
        'radius': 1.0,
        'ra_0': 0.0,
        'dec_0': 0.0, 
        'alpha': 30.0,
    }
}

RING = {
    'coords': 'ring',
    'params': {
        'seed_coord': 42,
        'radius': np.nan,
        'ra_0': 0.0,
        'dec_0': 0.0, 
        'alpha': 30.0,
    }
}