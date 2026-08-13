import numpy as np

F0_PSR, F1_PSR = 1.64e0, -4.69e-15
F0_MSC, F1_MSC = 2.72e2, -1.02e-15

def nan_freqs(n_array):
    F0 = np.full(n_array, fill_value=np.nan)
    F1 = np.full(n_array, fill_value=np.nan)
    return F0, F1

def constant_freqs(n_array, seed_freq = 42, f0=1.0, f1=np.nan):
    F0 = f0 * np.ones(n_array)
    F1 = f1 * np.ones(n_array)
    return F0, F1

FREQS = {
    'nan': {
        'distr': nan_freqs,
        'params': [],
    },
    'const': {
        'distr': constant_freqs,
        'params': ['seed_freq', 'f0', 'f1']
    }
}
