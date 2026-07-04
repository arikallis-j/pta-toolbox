import ptatoolbox as pta

# Create experiment directory
dm = pta.DataManager()
path = dm.create_experiment("test")
storage = dm.storage

# Load pre-cut atnf catalog 
cut_atnf = pta.load_catalog(storage, pta.CUT_ATNF_STEM, prefix=False)
print(cut_atnf.sample(n_psr=100, seed=42))

# Make synthetic mixed with pre-cut atnf catalog by method `cone`
catalog = pta.make_catalog(
    real_catalog=cut_atnf,
    n_psr=100,
    seed=42,
    name='test', 
    method='cone',
    params={'seed_psr': 42, 'alpha': 20.0, 'ra_0': 0.0, 'dec_0': 0.0},
    fields=['PSRJ', 'RAJD', 'DECJD']
)

print(catalog)
pta.plot_catalog(catalog, path)
pta.plot_pulsars(catalog, path)