import matplotlib.pyplot as plt
import numpy as np

def get_cartesian(ra, dec, px):
    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)
    px = np.asarray(px, dtype=float)
    ra = np.deg2rad(ra)
    dec = np.deg2rad(dec)
    D = np.where((px > 0) & np.isfinite(px), 1 / px, np.nan)
    
    x = D * np.cos(dec) * np.cos(ra)
    y = D * np.cos(dec) * np.sin(ra)
    z = D * np.sin(dec)
    return x, y, z

def plot_catalog(catalog, path, show=False, save=True):
    ra, dec = catalog.data['RAJD'], catalog.data['DECJD']
    ra = np.where(ra>=180, ra-360, ra)
    ra *= np.pi/180
    dec *= np.pi/180
    plt.figure(figsize=(8,4))
    plt.subplot(projection="mollweide")
    plt.title(f"{catalog.name} pulsar catalog")
    plt.grid(True)
    plt.plot(ra, dec, 'o', markersize=2)
    if save:
        plt.savefig(path / f"{catalog.name}_cat.png", dpi=1000)
    if show:
        plt.show()
    plt.close()

def plot_pulsars(catalog, path, show=False, save=True):
    ra, dec, px = catalog.data['RAJD'], catalog.data['DECJD'], catalog.data['PX']
    x, y, z = get_cartesian(ra, dec, px)
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x, y, z, c='blue', s=20, label='PTA')
    ax.scatter(0, 0, 0, c='red', s=100, label='Observer')

    ax.set_box_aspect([1, 1, 1])

    ax.set_xlabel("x, kpc")
    ax.set_ylabel("y, kpc")
    ax.set_zlabel("z, kpc")
    ax.set_title("Distribution of Pulsar Array")
    ax.legend()
    
    if save:
        plt.savefig(path / f"{catalog.name}_pta.png", dpi=1000)
    if show:
        plt.show()
    plt.close()