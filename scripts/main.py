import ptatoolbox as pta

dm = pta.DataManager(root_dir='./data')
path = dm.create_experiment("test")
storage = dm.storage

# Isotropic distribution on a sphere with 1 kpc radius
sphere = pta.make_catalog(
    n_psr = 1000,
    name='sphere', 
    method='sphere',
    params={'seed_psr': 42, 'radius': 1.0}
)
print(sphere)
pta.plot_catalog(sphere, path)
pta.plot_pulsars(sphere, path)