import matplotlib.pyplot as plt
import numpy as np
import ptatoolbox as pta
import typer

def main(exp: str = 'test-cat', flat: bool = False, n_psr: int = 10, distr: str = 'sphere'):
    dm = pta.utils.DataManager()
    sim_path = dm.make_dir(path=dm.raw, name=exp)
    res_path = dm.make_dir(path=dm.processed, name=exp)

    if distr == 'sphere':
        cat_config = pta.SPHERE
    elif distr == 'ball':
        cat_config = pta.BALL
    elif distr == 'cap':
        cat_config = pta.CAP
    elif distr == 'cone':
        cat_config = pta.CONE
    elif distr == 'ring':
        cat_config = pta.RING
    else:
        cat_config = pta.SPHERE

    psrcat = pta.make_catalog(n_psr, config = cat_config)
    print(psrcat)

    phi, theta, rho = pta.get_coords(psrcat.data)
    print(phi)
    print(theta)
    print(rho)
    
    if flat:
        fig, axis = plt.subplots(figsize=(10, 5))
    else:
        fig, axis = plt.subplots(subplot_kw={'projection': 'mollweide'}, figsize=(10, 5))

    sky = pta.plot_sky(
        phi, theta, axis, 
        title='Pulsar distribution', 
        marker='o',
        size=20,
        hours=True,
    )

    plt.show()

    return None

if __name__ == '__main__':
    typer.run(main)