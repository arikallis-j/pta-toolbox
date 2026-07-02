import ptatoolbox as pta

import pandas as pd

dm = pta.DataManager(root_dir='./data')
path = dm.create_experiment("test")
storage = path.parent / ".storage"
atnf = pta.load_catalog(storage, 'atnf')
catalog = pta.make_catalog(key='atnf', params={'base': atnf})
print(catalog.data)

tempo_format = {
    'RAJD': 'RAJ',
    'DECJD': 'DECJ',
}
   
def tempo_data(data):
    return data.rename(columns=tempo_format)
