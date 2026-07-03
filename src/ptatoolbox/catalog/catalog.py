import numpy as np
import pandas as pd

from typing import NamedTuple

from .io import load_data, dump_data
from .pulsar import make_pulsars, atnf_format
from .models import make_synthetics

tempo_format = {
    'RAJD': 'RAJ',
    'DECJD': 'DECJ',
}

class Catalog:
    def __init__(self, data, name='test'):
        self.name = name
        self.data = data[atnf_format.values()]

    def pulsars(self):
        return make_pulsars(self.data)

    def tempo(self):
        return self.data.rename(columns=tempo_format)

    def __repr__(self):
        return f"{self.name} catalog\n" + repr(self.data)

    def __str__(self):
        return f"{self.name} catalog\n" + str(self.data)

def make_catalog(n_psr, name='sample', method='test', params={}):
    data = make_synthetics(n_psr, method, params)
    return Catalog(data, name=name)

def load_catalog(path, name, prefix=True):
    filename = name
    if prefix:
        filename += "_cat"
    data = load_data(path, name=filename)
    return Catalog(data, name)

def save_catalog(catalog, path, prefix=True):
    filename = catalog.name
    if prefix:
        filename += "_cat"
    dump_data(catalog.data, path, name=filename)


    
