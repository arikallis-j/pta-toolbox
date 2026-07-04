import ptatoolbox as pta

# Create experiment directory
dm = pta.DataManager(root_dir='./data')
path = dm.create_experiment("test-pickle-io")

# Load data from the path
data = pta.load_data(dm.storage, stem='cut_atnf')
print(data)

# Dump data to the path
pta.dump_data(data, path, stem='data')