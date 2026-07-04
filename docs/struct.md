# structure

├── __init__.py
├── core/
│   ├── data_manager.py      # DataManager
│   ├── constants.py         # atnf_format, tempo_format
│   └── types.py             # Pulsar (может остаться в pulsar.py)
├── io/
│   ├── __init__.py
│   ├── + manager.py           # DataManager (можно вынести в core)
│   ├── + pickle_io.py         # load_data, dump_data
│   └── + atnf_io.py           # load_atnf, download_atnf (используют DataManager)
├── catalog/
│   ├── __init__.py
│   ├── catalog.py           # класс Catalog
│   ├── + pulsar.py            # Pulsar, make_pulsars, make_data
│   ├── factory.py           # make_catalog, from_real_and_synthetic, etc.
│   └── funcs.py             # вспомогательные функции (координаты и т.п.)
├── models/
│   ├── __init__.py
│   ├── generators.py        # simple_ring, simple_sphere и т.д.
│   └── synthetic.py         # make_synthetics, make_synthetics_coords_only
└── visuals/