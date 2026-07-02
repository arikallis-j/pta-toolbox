import shutil, pathlib

class DataManager:
    def __init__(self, root_dir=None):
        if root_dir is None :
            root_dir = "./data"
        
        self.root = pathlib.Path(root_dir).expanduser().resolve()
        self.storage = self.root / ".storage"
        self.locals = pathlib.Path(__file__).resolve().parent.parent / ".locals"

        self.root.mkdir(exist_ok=True)
        self.storage.mkdir(exist_ok=True)

        for path in self.locals.iterdir():
            if path.is_file():
                self.copy_from_locals(path.name)

    def create_experiment(self, name):
        self.experiment = self.root / name
        self.experiment.mkdir(exist_ok=True)
        return self.experiment

    def get_storage_path(self, filename=""):
        target = self.storage / filename
        return target

    def copy_from_locals(self, filename):
        shutil.copy(self.locals / filename, self.storage / filename)
        return self.storage / filename