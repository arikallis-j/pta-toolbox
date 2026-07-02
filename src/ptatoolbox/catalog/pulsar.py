import pandas as pd

from typing import NamedTuple

class Pulsar(NamedTuple):
    name: str
    ra: float # deg
    dec: float # deg
    f0: float = 1.0 # Hz
    f1: float = float('nan') # s^-2
    pmra: float = float('nan') # mas yr^-1
    pmdec: float = float('nan') # mas yr^-1
    px: float = float('nan') # mas
    dm: float = float('nan') # cm^-3 pc

atnf_format = {
    'name': 'PSRJ',
    'ra': 'RAJD',
    'dec': 'DECJD',
    'f0': 'F0',
    'f1': 'F1',
    'pmra': 'PMRA',
    'pmdec': 'PMDEC',
    'px': 'PX',
    'dm': 'DM',
}

def make_pulsar(params):
    return Pulsar(**params)

def make_pulsars(data):
    pulsars = []
    for row in data.itertuples():
        params = {attr: getattr(row, col) for attr, col in atnf_format.items()}
        pulsar = make_pulsar(params)
        pulsars.append(pulsar)
    return pulsars

def make_data(pulsars):
    data = pd.DataFrame(pulsars).rename(columns=atnf_format)
    return data

def make_test_array(n_psr):
    pulsars = []
    for k in range(n_psr):
        psr = Pulsar(name=f"J{k}", ra=k/n_psr * 360, dec=0.0)
        pulsars.append(psr)
    return pulsars