import matplotlib.pyplot as plt
import numpy as np
import ptatoolbox as pta
import typer

def main(    
    exp: str = 'test',
    name: str = 'psrcat',
    flat: bool = False,
    show: bool = False,
    save: bool = False
):
    dm = pta.utils.DataManager()
    sim_path = dm.make_dir(path=dm.raw, name=exp)
    res_path = dm.make_dir(path=dm.processed, name=exp)

    pulsars = pta.utils.load_data(path=sim_path, stem=name)
    phi, theta, rho = pta.get_coords(pulsars)
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
        fig.savefig(res_path / f'{name}.png', dpi=1000)
    if show:
        plt.show()
    return None

if __name__ == '__main__':
    typer.run(main)