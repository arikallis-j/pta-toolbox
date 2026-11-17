"""Manages file system paths and local data for ptatoolbox"""

import shutil
import yaml
import xarray as xr
from pathlib import Path
from typing import Optional, Union

def parse_yaml(path):
    with open(Path(path)) as f:
        return yaml.safe_load(f)

def save_dataset(data, path):
    data.to_netcdf(path, engine="h5netcdf")

def load_dataset(path):
    with xr.open_dataset(path, engine="h5netcdf") as data:
        return data.load()

class DataManager:
    """Manages file system paths and local data for ptatoolbox."""
    def __init__(self, root_dir: Optional[Union[str, Path]] = "."):
        """Initialize the DataManager."""
        self.root = Path(root_dir or ".").expanduser().resolve()
        self.data = self.root / "data"
        self.config = self.root / "config"
        self.raw = self.data / "raw"
        self.processed = self.data / "processed"

        self.root.mkdir(exist_ok=True)
        self.data.mkdir(exist_ok=True)
        self.config.mkdir(exist_ok=True)
        self.raw.mkdir(exist_ok=True)
        self.processed.mkdir(exist_ok=True)

    def make_dir(self, path, name: str) -> Path:
        """Create a new directory inside the given path."""
        new_dir = path / name
        new_dir.mkdir(exist_ok=True)
        return new_dir

    def make_path(self, path: str) -> Path:
        path = Path(path).expanduser().resolve()
        return path

    def make_experiment(self, name: str) -> Path:
        """Create a new experiment directory inside the root."""
        sim_path = self.make_dir(path=self.raw, name=name)
        res_path = self.make_dir(path=self.processed, name=name)
        return sim_path, res_path