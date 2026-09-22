"""Manages file system paths and local data for ptatoolbox"""

import shutil
from pathlib import Path
from typing import Optional, Union

class DataManager:
    """Manages file system paths and local data for ptatoolbox."""
    def __init__(self, root_dir: Optional[Union[str, Path]] = "./data"):
        """Initialize the DataManager."""
        self.root = Path(root_dir or "./data").expanduser().resolve()
        self.raw = self.root / "raw"
        self.processed = self.root / "processed"
        self.storage = self.processed / ".storage"

        self.root.mkdir(exist_ok=True)
        self.raw.mkdir(exist_ok=True)
        self.processed.mkdir(exist_ok=True)
        self.storage.mkdir(exist_ok=True)

        self._copy_raw_to_storage()

    def create_experiment(self, name: str) -> Path:
        """Create a new experiment directory inside the root."""
        experiment = self.processed / name
        experiment.mkdir(exist_ok=True)
        return experiment

    def make_dir(self, path, name: str) -> Path:
        """Create a new directory inside the given path."""
        new_dir = path / name
        new_dir.mkdir(exist_ok=True)
        return new_dir

    def storage_file_path(self, filename: str = "") -> Path:
        """Return a path to a file inside the storage directory."""
        return self.storage / filename

    def _copy_raw_to_storage(self) -> None:
        """Copy all regular files from locals directory to storage."""
        for path in self.raw.iterdir():
            if path.is_file():
                self._copy_from_raw(path.name)

    def _copy_from_raw(self, filename: str) -> Path:
        """Copy a single file from locals to storage."""
        src = self.raw / filename
        dst = self.storage / filename
        shutil.copy(src, dst)