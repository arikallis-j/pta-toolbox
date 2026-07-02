import numpy as np
import pandas as pd

from typing import NamedTuple
from .io import *
from .pulsar import *
from .synth_data import *


class Catalog:
    def __init__(self, data, name='test'):
        self.name = name
        self.data = data[atnf_format.values()]
        self.pulsars = make_pulsars(self.data)

def save_catalog(catalog, path, prefix=True):
    filename = catalog.name
    if prefix:
        filename += "_cat"
    dump_data(catalog.data, path, name=filename)

def load_catalog(path, name, prefix=True):
    filename = name
    if prefix:
        filename += "_cat"
    data = load_data(path, name=filename)
    return Catalog(data, name)

def make_catalog(name='sample', key='test', params={}):
    data = make_synthetics[key](**params)
    return Catalog(data, name=name)

    
