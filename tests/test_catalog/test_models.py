import ptatoolbox as pta

# Create experiment directory
dm = pta.DataManager(root_dir='./data')
path = dm.create_experiment("test-models")
storage = dm.storage

# Make test synthetic catalog
test = pta.make_catalog(
    n_psr = 10,
    name='simple', 
    method='test',
    params={}
)
print(test)

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

# Isotropic distribution in a ball with 1 kpc radius
ball = pta.make_catalog(
    n_psr = 1000,
    name='ball', 
    method='ball',
    params={'seed_psr': 42, 'radius': 1.0}
)
print(ball)
pta.plot_catalog(ball, path)
pta.plot_pulsars(ball, path)

# Isotropic distribution in a cone with 30.0 deg angle radius
cone = pta.make_catalog(
    n_psr = 1000,
    name='cone', 
    method='cone',
    params={'seed_psr': 42, 'radius': 1.0,  'alpha': 30.0, 'ra_0': 0.0, 'dec_0': 0.0}
)
print(cone)
pta.plot_catalog(cone, path)
pta.plot_pulsars(cone, path)

# Isotropic distribution on a cap with 30.0 deg angle radius
cap = pta.make_catalog(
    n_psr = 1000,
    name='cap', 
    method='cap',
    params={'seed_psr': 42, 'alpha': 30.0, 'ra_0': 0.0, 'dec_0': 0.0}
)
print(cap)
pta.plot_catalog(cap, path)

# Isotropic distribution on a ring with 30.0 deg angle radius
ring = pta.make_catalog(
    n_psr = 1000,
    name='ring', 
    method='ring',
    params={'seed_psr': 42, 'alpha': 30.0, 'ra_0': 0.0, 'dec_0': 0.0}
)
print(ring)
pta.plot_catalog(ring, path)