"""Manages file system paths and local data for ptatoolbox"""

import shutil
from pathlib import Path
from typing import Optional, Union

class DataManager:
    """Manages file system paths and local data for ptatoolbox."""
    def __init__(self, root_dir: Optional[Union[str, Path]] = None):
        """Initialize the DataManager."""
        self.root = Path(root_dir or "./data").expanduser().resolve()
        self.storage = self.root / ".storage"
        self.root.mkdir(exist_ok=True)
        self.storage.mkdir(exist_ok=True)

        locals_dir = Path(__file__).resolve().parent.parent / ".locals"
        self.locals = Path(locals_dir).expanduser().resolve()
        self._copy_locals_to_storage()

    def create_experiment(self, name: str) -> Path:
        """Create a new experiment directory inside the root."""
        experiment = self.root / name
        experiment.mkdir(exist_ok=True)
        return experiment

    def storage_file_path(self, filename: str = "") -> Path:
        """Return a path to a file inside the storage directory."""
        return self.storage / filename

    def _copy_locals_to_storage(self) -> None:
        """Copy all regular files from locals directory to storage."""
        for path in self.locals.iterdir():
            if path.is_file():
                self._copy_from_locals(path.name)

    def _copy_from_locals(self, filename: str) -> Path:
        """Copy a single file from locals to storage."""
        src = self.locals / filename
        dst = self.storage / filename
        shutil.copy(src, dst)