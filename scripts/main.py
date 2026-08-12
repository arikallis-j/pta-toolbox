import ptatoolbox as pta

# Create experiment directory
dm = pta.DataManager()
path = dm.create_experiment("test")
storage = dm.storage

atnf = pta.load_catalog(storage, pta.ATNF_STEM, prefix=False)
cut_atnf = pta.load_catalog(storage, pta.CUT_ATNF_STEM, prefix=False)
cut_atnf.name = 'atnf'

config = {
    'coords': 'sphere',
    # 'freqs': 'const',
    # 'params': {
    #     'f0': pta.F0_PSR,
    #     'f1': pta.F1_PSR,
    # }
}

simple_cat = pta.make_catalog(100, config, mode='simple')
print(simple_cat)
pta.plot_catalog(simple_cat, path)