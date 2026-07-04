"""Shared constants and mappings for ptatoolbox."""

# Mapping from Pulsar attribute names to ATNF column names
ATNF_FORMAT = {
    'name': 'PSRJ',
    'ra': 'RAJD',
    'dec': 'DECJD',
    'f0': 'F0',
    'f1': 'F1',
    'pmra': 'PMRA',
    'pmdec': 'PMDEC',
    'px': 'PX',
    'dm': 'DM',
}

# Mapping from ATNF column names to Tempo2 column names
TEMPO_FORMAT = {
    'RAJD': 'RAJ',
    'DECJD': 'DECJ',
}
