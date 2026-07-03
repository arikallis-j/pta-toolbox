import psrqpy
import pandas as pd

# Data processing

def load_data(path, name):
    data = pd.read_pickle(path / f"{name}.pkl")
    return data

def dump_data(data, path, name='data'):
    data.to_pickle(path / f"{name}.pkl")

# ATNF working

def load_atnf(path):
    storage = path.parent / ".storage"
    atnf = load_data(storage, 'atnf')
    return atnf

def download_atnf(path):
    atnf = psrqpy.QueryATNF().pandas
    dump_data(atnf, path, 'atnf')
    return atnf

def load_cut_atnf(path):
    storage = path.parent / ".storage"
    atnf_cut = load_data(storage, 'cut_atnf')
    return atnf_cut

def cut_atnf(atnf, template):
    atnf_cut = atnf[atnf['PSRJ'].isin(template['PSRJ'])].reset_index()[template.columns]
    return atnf_cut
