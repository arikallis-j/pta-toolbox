import ptatoolbox as pta 
import pta_catalog as cat

# Create experiment directory
dm = pta.DataManager()
path = dm.create_experiment("test")
storage = dm.storage

atnf = cat.load_catalog(storage, cat.ATNF_STEM, prefix=False)
cut_atnf = cat.load_catalog(storage, cat.CUT_ATNF_STEM, prefix=False)
cut_atnf.name = 'atnf'

config = {
    'coords': 'sphere',
}

simple_cat = cat.make_catalog(100, config, mode='simple')
print(simple_cat)
cat.plot_catalog(simple_cat, path)