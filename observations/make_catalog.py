import ptatoolbox as pta
import pta_catalog as ptacat
import numpy as np
import typer
import os, glob

from enterprise.pulsar import Pulsar, Tempo2Pulsar

PARAMS = {
    'RAJ': 'ra',
    'DECJ': 'dec',
    'F0': 'f0',
    'F1': 'f1',
    'PMRA': 'pmra',
    'PMDEC': 'pmdec',
    'PX': 'px',
    'DM': 'dm',
}

def main(nside: int = 100, sample = "psr_test_1"):
    lmax = 2 * nside

    dm = pta.utils.DataManager()
    path = dm.create_experiment("observations")
    datapath = dm.raw / sample
    parpath = datapath / "par" / "*.par"
    timpath = datapath / "tim" / "*.tim"

    parfiles = sorted(glob.glob(str(parpath)))
    timfiles = sorted(glob.glob(str(timpath)))

    PSRs = []
    pulsars = []

    for p, t in zip(parfiles, timfiles):
        psr = Pulsar(p, t, timing_package='tempo2', drop_t2pulsar=False)
        PSRs.append(psr)
        
        params = {PARAMS[key]: psr.t2pulsar[key].val for key in PARAMS.keys()}
        params['origin'] = 'SIM'
        params['name'] = psr.name
        params['ra'] = np.rad2deg(params['ra'])
        params['dec'] = np.rad2deg(params['dec'])
        pulsar = ptacat.make_pulsar(params)
        pulsars.append(pulsar)

    data = ptacat.make_data(pulsars)
    psrcat = ptacat.Catalog(data=data, name=sample)
    ptacat.save_catalog(psrcat, datapath)
    print(psrcat)    

if __name__ == '__main__':
    typer.run(main)