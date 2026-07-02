from .io import *
from .pulsar import *

def make_test(n_psr):
    pulsars = []
    for k in range(n_psr):
        psr = Pulsar(name=f"J{k}", ra=k/n_psr * 360, dec=0.0)
        pulsars.append(psr)
    data = make_data(pulsars)
    return data

def make_atnf(base):
    data = base.data
    return data

make_synthetics = {
    'test': make_test,
    'atnf': make_atnf,
}