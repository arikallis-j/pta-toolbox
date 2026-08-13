# Синтетический пульсарный каталог

## Формирование сигнала

1. Частотные характеристики: f0, f1
    > Определяют частотный портрет сигнала
    - Synthetic Periods: одинаковые f0 и f1, LogUniform f0, Normal f0, Binormal f0 (два класса пульсаров - обычные и милисекундные, например), Spacial Correlations (тест на ложную анизотропию)
    - Period Models: normal, log-normal, MSP model, disk/slab/isotropic, radial profile, scale height

    Synthetic Periods - не реализован, нужна кастомная версия
    Period Models - не реализован, практически всё есть в библиотеке `PsrPopPy`

2. Координаты: ra, dec
   > Определяют геометрию сигнала
   - Synthetic geometry: сфера, шапка, кольцо
   - Galactic population: Lorimer, Yusifov & Küçük, Faucher-Giguère & Kaspi, Spiral arm models, Kick velocity

   Synthetic geometry - реализована кастомная версия
   Galactic population - не реализован, практически всё есть в библиотеке `PsrPopPy`

## Астрометрия

1. Параллакс: px
   > Влияет на dm и на f1 через Shklovskii effect
   - Synthetic parallax: без параллакса, одинаковый px (сфера, шапка, кольцо)
   - Galactic parallax: Lorimer, Yusifov & Küçük, Faucher-Giguère & Kaspi, Spiral arm models, Kick velocity

   Synthetic parallax - реализована кастомная версия
   Galactic parallax - не реализован, практически всё есть в библиотеке `PsrPopPy`

2. Собственное движение: pmra, pmdec
    > Влияет на f1 через Shklovskii effect
    - Synthetic Proper Motion: без собственного движения, Spacial Correlations (тест на ложную анизотропию)
    - Galactic Proper Motion: kick velocity (Maxwellian), velocity model

    Synthetic Proper Motion - не реализовано, нужна кастомная версия
    Galactic Proper Motion - не реализован, практически всё есть в библиотеке `PsrPopPy`

## Искажение сигнала

1. Мера дисперсии dm
    > Дает частотно зависимые задержки
    - Synthetic Dispersion: без меры дисперсии, Spacial Correlations (тест на ложную анизотропию)
    - Galactic Dispersion: YMW16, NE2001, NE2025

    Synthetic Dispersion - не реализовано, нужна кастомная версия
    Galactic Dispersion - не реализован, практически всё есть в библиотеке `PyGEDM`

`PsrPopPy`

- git: https://github.com/samb8s/PsrPopPy
- doc: https://samb8s.github.io/PsrPopPy/

```bash
conda activate pta
pip install enscons
python setup.py install
# Some errors
pip install --prefix=$CONDA_PREFIX dist/psrpoppy-1.0.0-cp311-cp311-linux_x86_64.whl
```

`PyGEDM`

- git: https://github.com/FRBs/pygedm
- doc: https://pygedm.readthedocs.io/en/latest/
