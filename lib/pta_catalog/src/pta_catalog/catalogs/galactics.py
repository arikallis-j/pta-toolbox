import numpy as np
import pandas as pd
from .coords import COORDS, isotropic_ball
from .freqs import FREQS, F0_PSR, F1_PSR

def get_params(config, params):
    return {
        key: config[key] 
        for key in config.keys()
        if key in params
    }

def simple_galactic(n_array: int, config = {}):
    coord = COORDS[config.get('coords', 'sphere')]
    freqs = FREQS[config.get('freqs', 'nan')]
    params = config.get('params', {})
    Ra, Dec, Px = coord['distr'](n_array, **get_params(params, coord['params']))
    F0, F1  = freqs['distr'](n_array, **get_params(params, freqs['params']))
    df = pd.DataFrame(np.array([Ra, Dec, Px, F0, F1]).T, columns=['RAJD','DECJD', 'PX', 'F0', 'F1'])
    df['TYPE'] = 'SIM'
    return df

def physical_galactic(n_array: int, config = {}):
    coord = COORDS['ball']
    freqs = FREQS['const']
    params = config.get('params', {})
    Ra, Dec, Px = coord['distr'](n_array, **get_params(params, coord['params']))
    F0, F1  = freqs['distr'](n_array, **get_params(params, freqs['params']))
    df = pd.DataFrame(np.array([Ra, Dec, Px, F0, F1]).T, columns=['RAJD','DECJD', 'PX', 'F0', 'F1'])
    df['TYPE'] = 'SYN'
    return df

GALACTICS = {
    'simple': simple_galactic,
    'physical': physical_galactic,
}