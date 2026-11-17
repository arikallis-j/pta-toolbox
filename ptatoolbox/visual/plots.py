import numpy as np
import healpy as hp
import matplotlib
import matplotlib.pyplot as plt

# from matplotlib.colors import LogNorm

DEG = np.pi/180

def map_decomposition(h):
    amplitude = np.abs(h)
    phase = np.mod(np.angle(h), 2 * np.pi)
    return amplitude, phase

def plot_sky(phi, theta, axis, title='Sky', color='steelblue', marker='o', grid=True, size=10, alpha=0.7, hours=False):
    lon = 2*np.pi - phi
    lon = np.mod(lon + np.pi, 2 * np.pi) - np.pi
    lat = np.pi/2.0 - theta
    sky = axis.scatter(lon, lat, s=size, alpha=alpha, color=color, marker=marker)
    axis.set_title(title)
    if grid:
        axis.grid(True, alpha=0.3)
    axis.set_xlabel("RA")
    axis.set_ylabel("Dec")
    axis.set_yticks([-75*DEG, -60*DEG, -45*DEG, -30*DEG, -15*DEG, 0*DEG, +15*DEG, +30*DEG, +45*DEG, +60*DEG, +75*DEG], ['-75°', '-60°', '-45°', '-30°', '-15°', '0°', '+15°', '+30°', '+45°', '+60°', '+75°'])
    axis.set_xticks([-150*DEG, -120*DEG, -90*DEG, -60*DEG, -30*DEG, 0*DEG, +30*DEG, +60*DEG, +90*DEG, +120*DEG, +150*DEG], ['150°', '120°', '90°', '60°', '30°', '0°', '330°', '300°', '270°', '240°', '210°'])
    if hours:
        axis.set_xticks([-150*DEG, -120*DEG, -90*DEG, -60*DEG, -30*DEG, 0*DEG, +30*DEG, +60*DEG, +90*DEG, +120*DEG, +150*DEG], ['10ʰ', '8ʰ', '6ʰ', '4ʰ', '2ʰ', '0ʰ', '22ʰ', '20ʰ', '18ʰ', '16ʰ', '14ʰ'])
    if not isinstance(axis, matplotlib.projections.geo.MollweideAxes):
        axis.set_xlim(-np.pi, np.pi)
        axis.set_ylim(-np.pi/2, np.pi/2)
    return sky

def plot_sky_healpix(map_pix, axis,  title='Sky', cmap='viridis', grid=True, alpha=1.0, vmin=None, vmax=None, n_grid=1000, hours=False):
    nside = hp.npix2nside(len(map_pix))
    lon = np.linspace(-np.pi, np.pi, 2*n_grid)
    lat = np.linspace(-np.pi / 2.0, np.pi / 2.0, n_grid)
    LON, LAT = np.meshgrid(lon, lat)

    theta = np.pi / 2.0 - LAT
    phi = np.mod(2 * np.pi - LON, 2 * np.pi)

    pix = hp.ang2pix(nside, theta, phi, nest=False)
    grid_data = map_pix[pix]

    sky = axis.pcolormesh(LON, LAT, grid_data, cmap=cmap, alpha=alpha, shading='auto', vmin=vmin, vmax=vmax, rasterized=True)
    axis.set_title(title)
    if grid:
        axis.grid(True, alpha=0.3)
    axis.set_xlabel("RA")
    axis.set_ylabel("Dec")
    axis.set_yticks([-75*DEG, -60*DEG, -45*DEG, -30*DEG, -15*DEG, 0*DEG, +15*DEG, +30*DEG, +45*DEG, +60*DEG, +75*DEG], ['-75°', '-60°', '-45°', '-30°', '-15°', '0°', '+15°', '+30°', '+45°', '+60°', '+75°'])
    axis.set_xticks([-150*DEG, -120*DEG, -90*DEG, -60*DEG, -30*DEG, 0*DEG, +30*DEG, +60*DEG, +90*DEG, +120*DEG, +150*DEG], ['150°', '120°', '90°', '60°', '30°', '0°', '330°', '300°', '270°', '240°', '210°'])
    if hours:
        axis.set_xticks([-150*DEG, -120*DEG, -90*DEG, -60*DEG, -30*DEG, 0*DEG, +30*DEG, +60*DEG, +90*DEG, +120*DEG, +150*DEG], ['10ʰ', '8ʰ', '6ʰ', '4ʰ', '2ʰ', '0ʰ', '22ʰ', '20ʰ', '18ʰ', '16ʰ', '14ʰ'])
    if not isinstance(axis, matplotlib.projections.geo.MollweideAxes):
        axis.set_xlim(-np.pi, np.pi)
        axis.set_ylim(-np.pi/2, np.pi/2)
    return sky

def plot_sky_amplitude(map_pix, axis, title='Amplitude', cmap='viridis', grid=True, alpha=1.0, hours=False, vmax=None, n_grid=1000):
    amplitude = np.abs(map_pix)
    sky = plot_sky_healpix(amplitude, axis, title=title, cmap=cmap, grid=grid, alpha=alpha, hours=hours, vmin=0.0, vmax=vmax, n_grid=n_grid)
    return sky

def plot_sky_phase(map_pix, axis, title='Phase', cmap='hsv', grid=True, alpha=1.0, hours=False, n_grid=1000):
    phase = np.mod(np.angle(map_pix), 2 * np.pi)
    sky = plot_sky_healpix(phase, axis, title=title, cmap=cmap, grid=grid, alpha=alpha, hours=hours, vmin=0.0, vmax=2*np.pi, n_grid=n_grid)
    return sky

def plot_sky_real(map_pix, axis, title='Real', cmap='plasma', grid=True, alpha=1.0, hours=False, vmin=None, vmax=None, n_grid=1000):
    real = np.real(map_pix)
    if vmin is None and vmax is None:
        vabs = np.max(np.abs(real))
        vmin, vmax = -vabs, +vabs
    sky = plot_sky_healpix(real, axis, title=title, cmap=cmap, grid=grid, alpha=alpha, hours=hours, vmin=vmin, vmax=vmax, n_grid=n_grid)
    return sky

def plot_sky_imag(map_pix, axis, title='Imag', cmap='plasma', grid=True, alpha=1.0, hours=False, vmin=None, vmax=None, n_grid=1000):
    imag = np.imag(map_pix)
    if vmin is None and vmax is None:
        vabs = np.max(np.abs(imag))
        vmin, vmax = -vabs, +vabs
    sky = plot_sky_healpix(imag, axis, title=title, cmap=cmap, grid=grid, alpha=alpha, hours=hours, vmin=vmin, vmax=vmax, n_grid=n_grid)
    return sky
    
def plot_colorbar_healpix(sky, fig, axis, orientation='horizontal', shrink=0.5, aspect=30, pad=0.1):
    cbar = fig.colorbar(sky, ax=axis, orientation=orientation, shrink=shrink, aspect=aspect, pad=pad)
    return cbar


PLOTS = {
    'sky': plot_sky,
}

PLOTS_COLORBAR = {
    'sky': True,
}
