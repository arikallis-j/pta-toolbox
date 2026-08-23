import ptatoolbox as pta
import numpy as np
import typer

def main(nside: int = 50, polar: str = 'n'):
    dm = pta.utils.DataManager()
    path = dm.create_experiment("stokes")

    alm_left = pta.generate_flm(nside, key=42, L_lower=2)
    alm_right = pta.generate_flm(nside, key=43, L_lower=2)

    h_left = pta.inverse_sfft2(alm_left, nside, spin=2)
    h_right = pta.inverse_sfft2(alm_right, nside, spin=-2)
    h_plus, h_cross = pta.circ2lin(h_left, h_right) 

    if polar == 'l':   # circular 'l' polarization
        h_plus, h_cross = pta.circ2lin(h_left, np.zeros_like(h_right)) 
    elif polar == 'r': # circular 'r' polarization
        h_plus, h_cross = pta.circ2lin(np.zeros_like(h_left), h_right) 
    elif polar == '+': # linear '+'
        h_plus, h_cross = h_plus, np.zeros_like(h_cross) 
    elif polar == 'x': # linear 'x'
        h_plus, h_cross = np.zeros_like(h_plus), h_cross
    else:              # no polarization
        h_plus, h_cross = h_plus, h_cross

    I, Q, U, V = pta.Stokes_parameters(h_plus, h_cross)
    I_l, I_c = pta.polarizarion_intensities(h_plus, h_cross)
    psi, chi = pta.polarizarion_angles(h_plus, h_cross)

    pta.plot_map(h_plus, name='h-plus', path=path, show=False)
    pta.plot_map(h_cross, name='h-cross', path=path, show=False)
    
    pta.plot_intensities(h_plus, h_cross, path=path, show=False)
    pta.plot_polarization_degrees(h_plus, h_cross, path=path, show=False)
    pta.plot_Stokes_parameters(h_plus, h_cross, path=path, show=False)

    P_deg = pta.polarizarion_degree(h_plus, h_cross)
    print(f"Polarizarion degree: {P_deg:.2e}")

if __name__ == '__main__':
    typer.run(main)