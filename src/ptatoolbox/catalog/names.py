import numpy as np
from collections import Counter

def get_name(ra, dec, prefix='S'):
    ra_h = np.floor(ra/15.0).astype(int)
    ra_m = np.floor(np.round((ra/15.0 - ra_h) * 60, decimals=1)).astype(int)
    dec_d = np.floor(np.abs(dec)).astype(int)
    dec_m = np.floor(np.round((np.abs(dec) - dec_d) * 60, decimals=1)).astype(int)
    dec_s = "+" if np.sign(dec)>=0.0 else "-"
    name = f"{prefix}{ra_h:02d}{ra_m:02d}{dec_s}{dec_d:02d}{dec_m:02d}"
    return name

def get_names(ra, dec, prefix='S'):
    base_names = [get_name(r, d, prefix) for r, d in zip(ra, dec)]
    freq = Counter(base_names)
    counters = {name: 0 for name in freq}
    final_names = []
    for name in base_names:
        cnt = counters[name]
        if cnt == 0:
            final_names.append(name)
        else:
            if cnt <= 26:
                suffix = chr(ord('A') + cnt - 1)
            else:
                idx = cnt - 27
                first = chr(ord('a') + idx // 26)
                second = chr(ord('a') + idx % 26)
                suffix = first + second
            final_names.append(name + suffix)
        counters[name] += 1
    
    return final_names