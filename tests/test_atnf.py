import ptatoolbox as pta

# Create experiment directory
dm = pta.DataManager(root_dir='./data')
path = dm.create_experiment("test-atnf")
storage = dm.storage

# Load the full ATNF DataFrame from storage
atnf = pta.load_atnf(path)
print(atnf)

# Load cut ATNF DataFrame from storage
atnf_cut = pta.load_cut_atnf(path)
print(atnf_cut)

# Download the full ATNF DataFrame by `psrqpy`
new_atnf = pta.download_atnf(path)
print(new_atnf)

# Cut ATNF DataFrame by a template
new_atnf_cut = pta.cut_atnf(new_atnf, atnf_cut)
print(new_atnf_cut)