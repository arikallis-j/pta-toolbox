import ptatoolbox as pta

# Create root directory in `data`
dm = pta.DataManager(root_dir='./data')
print(dm.root)
print(dm.storage)
print(dm.locals)

# Find files in storage
atnf_path = dm.get_storage_path(filename='atnf_cat.pkl')
print(atnf_path)

# Create experiment directory in `test-io-manager`
path = dm.create_experiment("test-io-manager")
print(path)