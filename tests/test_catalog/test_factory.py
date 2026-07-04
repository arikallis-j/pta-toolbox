import ptatoolbox as pta

# Create experiment directory
dm = pta.DataManager(root_dir='./data')
path = dm.create_experiment("test-factory")
storage = dm.storage

# Make full synthetic catalog by method `test`
catalog = pta.make_catalog(
    n_psr = 10,
    name='simple', 
    method='test',
    params={}
)
print(catalog)

# Save new catalog with prefix `_cat`
pta.save_catalog(catalog, path, prefix=True)

# Load pre-cut atnf catalog 
cut_atnf = pta.load_catalog(storage, pta.CUT_ATNF_STEM, prefix=False)
print(cut_atnf)

# Sample pre-cut atnf catalog
sample = cut_atnf.sample(n_psr=10, seed=42)
print(sample)

# Make synthetic mixed with real catalog by method `test`
catalog = pta.make_catalog(
    n_psr = 10,
    real_catalog=cut_atnf,
    name='new',
    method='test',
    params={},
    seed=42,
    fields=['PSRJ', 'RAJD', 'DECJD']
)
print(catalog)

