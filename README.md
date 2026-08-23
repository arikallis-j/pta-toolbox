# pta-toolbox

Pulsar Timing Array Toolbox for GW simulations and sky imaging

## Setup enviroment

Create conda enviroment:

```bash
conda env create -f environment/conda_env.yaml
```

Update conda enviroment:

```bash
conda env update -f environment/conda_env.yaml
```

Download pip dependencies:

```bash
cd environment && pip install -r requirements.txt && cd ..
```

Freeze current python dependencies:

```bash
pip list --not-required --format=freeze > environment/backup/requirements.txt.backup
```
