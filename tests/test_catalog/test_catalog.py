import ptatoolbox as pta

# Create experiment directory
dm = pta.DataManager(root_dir='./data')
path = dm.create_experiment("test-catalog")
storage = dm.storage

# Make empty catalog
catalog = pta.Catalog(name='small')
print(catalog)

# Add pulsars into new catalog
psr_1 = pta.Pulsar(name=pta.get_name(90.0, 0.0), ra=90.0, dec=0.0)
psr_2 = pta.Pulsar(name=pta.get_name(270.0, 0.0), ra=270.0, dec=0.0)
catalog = catalog.add([psr_1, psr_2])
print(catalog)

# Catalog's parameters and methods
print(catalog.name)
print(catalog.data)
print(catalog.tempo())
print(catalog.pulsars())

# Make new catalog by the data
data = pta.load_data(storage, 'cut_atnf')
catalog = pta.Catalog(data=data, name='cut_atnf')
print(catalog)

# Plot pulsar's map for catalog
pta.plot_catalog(catalog, path)