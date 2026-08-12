import numpy as np
import astropy.units as u

N = 10000

# радиус рождения в диске
# экспоненциальный диск
R0 = 4.5 # kpc

R = np.random.gamma(
    shape=2,
    scale=R0,
    size=N
)

# ограничим разумным радиусом
R = np.clip(R, 0.1, 20)


# азимут
phi = np.random.uniform(
    0,
    2*np.pi,
    N
)


# высота над диском
z = np.random.normal(
    0,
    0.05,
    N
) # kpc


sigma = 265 # km/s

vx_kick = np.random.normal(
    0,
    sigma,
    N
)

vy_kick = np.random.normal(
    0,
    sigma,
    N
)

vz_kick = np.random.normal(
    0,
    sigma,
    N
)

from galpy.orbit import Orbit
from galpy.potential import MWPotential2014


vo = 232
ro = 8.2


# скорости вращения диска
vR = vx_kick
vT = 220 + vy_kick
vz_vel = vz_kick


orb = Orbit(
    [
        R,
        vR,
        vT,
        z,
        vz_vel,
        phi
    ],
    ro=ro,
    vo=vo
)