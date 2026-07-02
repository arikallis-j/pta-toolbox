import psrqpy
import pandas as pd

def load_data(path, name):
    data = pd.read_pickle(path / f"{name}.pkl")
    return data

def dump_data(data, path, name='data'):
    data.to_pickle(path / f"{name}.pkl")

def load_atnf(path, update=False):
    if update:
        update_atnf(path)
    storage = path.parent / ".storage"
    atnf = pd.read_pickle(storage / "atnf.pkl")
    return atnf

def load_atnf_catalog(path, update=False):
    if update:
        update_atnf_catalog(path)
    storage = path.parent / ".storage"
    catalog = pd.read_pickle(storage / "atnf_cat.pkl")
    return catalog

def update_atnf(path):
    storage = path.parent / ".storage"
    atnf = psrqpy.QueryATNF().pandas
    atnf.to_pickle(storage / "atnf.pkl")

def update_atnf_catalog(path):
    storage = path.parent / ".storage"
    catalog = load_atnf_catalog(path)
    atnf = update_atnf(path)
    catalog = atnf[atnf['PSRJ'].isin(catalog['PSRJ'])].reset_index()[catalog.columns]
    catalog.to_pickle(storage / "atnf_cat.pkl")
