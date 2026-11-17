import matplotlib.pyplot as plt
import numpy as np
import ptatoolbox as pta
import typer

def main(
    exp: str = 'test',
    name: str = 'psrcat',
    n_psr: int = 10, 
    distr: str = 'sphere',
    seed: int = 42,
    radius: float = np.nan,
    ra_0: float = 0.0,
    dec_0: float = 0.0,
    alpha: float = 30.0,
    flat: bool = False,
    show: bool = False,
    save: bool = False,
):
    dm = pta.utils.DataManager()
    sim_path = dm.make_dir(path=dm.raw, name=exp)
    res_path = dm.make_dir(path=dm.processed, name=exp)

    config = {
        'coords': distr,
        'params': {
            'seed_coord': seed,
            'radius': radius,
            'ra_0': ra_0,
            'dec_0': dec_0,
            'alpha': alpha,
        }
    }

    for key, value in config.items():
        print(f"{key}:")
        if isinstance(value, dict):
            for k, v in value.items():
                print(f"\t{k:10} {v}")
        else:
            print(f"\t{value}")
    pulsars = pta.make_pulsars(n_psr, config=config)
    pta.utils.dump_data(data=pulsars, path=sim_path, stem=name)
    print(pulsars)

    phi, theta, rho = pta.get_coords(pulsars)
    
    if flat:
        fig, axis = plt.subplots(figsize=(10, 5))
    else:
        fig, axis = plt.subplots(subplot_kw={'projection': 'mollweide'}, figsize=(10, 5))

    sky = pta.plot_sky(
        phi, theta, axis, 
        title='Pulsar distribution', 
        marker='o',
        size=10,
        hours=True,
    )

    if save:
        fig.savefig(sim_path / f'{name}.png', dpi=1000)
    if show:
        plt.show()

    return None

if __name__ == '__main__':
    typer.run(main)