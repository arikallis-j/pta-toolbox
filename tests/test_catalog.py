import ptatoolbox as pta

# Create experiment directory
dm = pta.DataManager(root_dir='./data')
path = dm.create_experiment("test-catalog")
storage = dm.storage

# Make simple catalog by method `test`
catalog = pta.make_catalog(
    n_psr = 10,
    name='simple', 
    method='test',
    params={}
)
print(catalog)

# Catalog's parameters and methods
print(catalog.name)
print(catalog.data)
print(catalog.tempo())
print(catalog.pulsars())

# Change catalog's data
catalog.data = catalog.data[::2]
print(catalog)

# Save new catalog with prefix `_cat`
pta.save_catalog(catalog, path, prefix=True)

# Load new catalog
catalog = pta.load_catalog(path, 'simple', prefix=True)
print(catalog)

# Load full ATNF Catalog from storage
atnf = pta.load_catalog(storage, 'atnf', prefix=False)
print(atnf)

# Make new catalog by the data
data = pta.load_data(storage, 'cut_atnf')
cut_atnf = pta.Catalog(data=data, name='cut_atnf')
print(cut_atnf)

# Plot pulsar's map for catalog
pta.plot_catalog(atnf, path)
