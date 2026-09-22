import ptatoolbox as pta
# import matplotlib.pyplot as plt
# import numpy as np
import xarray as xr
import typer

def main(n_freqs: int = -1, sample: str = "psr_test_1"):
    dm = pta.utils.DataManager()
    path = dm.create_experiment("observations")
    datapath = dm.raw / sample

    pulsars = pta.utils.load_pulsars(path=datapath)
    tspan, freqs, n_freqs = pta.utils.global_frequency_grid(pulsars, n_freqs=n_freqs)
    print(f"tspan[years]: {tspan/(3600*24*365.25):.2f}")
    print(f"n_freqs: {n_freqs}")
    print(f"freqs[nHz]: {freqs/1e-9}")

    freqs, r_freqs, N_covs = pta.utils.project_to_fourier_domain(pulsars, n_freqs=n_freqs)
    print(f"freqs.shape: {freqs.shape}")
    print(f"r_freqs.shape: {r_freqs.shape}")
    print(f"N_covs.shape: {N_covs.shape}")

    theta, phi = pta.utils.load_pulsar_coords(pulsars)
    names = pta.utils.load_pulsar_names(pulsars)
    print(f"thetas.shape: {theta.shape}")
    print(f"phis.shape: {phi.shape}")
    print(f"names.shape: {names.shape}")

    ds = xr.Dataset(
        {
            "r_freq": (["freq", "psr"], r_freqs),
            "N_cov":  (["freq", "psr", "psr"], N_covs),
        },
        coords={
            "freq":   freqs,
            "psr":    names,
            "theta": (["psr"], theta),  
            "phi":   (["psr"], phi),   
        },
        attrs={
            "tspan":   tspan,
            "n_freqs": n_freqs,
            "sample":  sample,
        }
    )
    
    ds.to_netcdf(path / f"{sample}.nc", engine="h5netcdf")

    with xr.open_dataset(path / f"{sample}.nc" , engine="h5netcdf") as data:
        ds = data.load()

    print(ds.N_cov.shape)


if __name__ == '__main__':
    typer.run(main)

    # psrcat = ptacat.load_catalog(datapath, name=sample)
    # print(psrcat)

    # psr_name = "J100"
    # freq = ds.freq
    # r_freq = ds.r_freq.sel(pulsar=psr_name)
    # spectrum = np.real(r_freq * np.conj(r_freq))
    # plt.loglog(freq/1e-9, spectrum)
    # plt.show()


    # idx_psr = 0
    # spectra = np.real(r_freq * np.conj(r_freq))
    # plt.loglog(freqs/1e-9, spectra[:, idx_psr])
    # plt.show()

    # test_PSR = pulsars[3]
    # print(test_PSR.name)
    # print(test_PSR.phi * 180/np.pi)
    # print((np.pi/2 - test_PSR.theta) * 180/np.pi)
    # plt.plot(test_PSR.toas/(24*3600), test_PSR.residuals)
    # plt.show()
