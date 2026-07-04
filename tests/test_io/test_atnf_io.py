import ptatoolbox as pta

# Create experiment directory
dm = pta.DataManager()
path = dm.create_experiment("test-atnf-io")
storage = dm.storage

# Download the full ATNF DataFrame by `psrqpy`
new_atnf = pta.download_atnf(path)
print(new_atnf)

# Load the full ATNF DataFrame from storage
atnf = pta.load_data(storage, pta.ATNF_STEM)
print(atnf)

# Load pre-cut ATNF DataFrame from storage
atnf_cut = pta.load_data(storage, pta.CUT_ATNF_STEM)
print(atnf_cut)

# Get DataFrame with only required columns from ATNF_FORMAT
data = atnf.drop(['PSRJ', 'DM'], axis=1)
normalized_atnf = pta.normalize_atnf_data(data)
print(normalized_atnf)