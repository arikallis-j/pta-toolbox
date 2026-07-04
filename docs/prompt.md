# prompt

Привет! У меня возник вопрос по поводу правильной архитектуры кода. Я пишу библиотеку `pta-toolbox`, в которой я хочу собрать все необходимые инструменты для работы с пульсарным таймингом и поиском гравитационных волн. На данный момент, структура проекта выглядит примерно так:
- src/ptatoolbox/:
- - .locals/
- - io/
- - catalog/
- - visuals/

Где в .locals у меня просто лежат необходимые файлы (небольшие) с данными, которые могут понадобится для работы, а в visuals просто лежат несколько функций для красивой отрисовки каких-то объектов, так что пока про них не будем говорить. Основаную часть кода занимают io и catalog. Для начала io. В этой папке у меня по сути лежит только один файл manager.py с вот таким классом:
```python
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
                self._copy_from_locals(path.name)

    def create_experiment(self, name):
        experiment = self.root / name
        experiment.mkdir(exist_ok=True)
        return experiment

    def get_storage_path(self, filename=""):
        target = self.storage / filename
        return target

    def _copy_from_locals(self, filename):
        shutil.copy(self.locals / filename, self.storage / filename)
        return self.storage / filename
```
Этот DataManager просто аккуратно работает с путями и создает папки для работы. Возможно сюда я добавлю одну функцию для создания папок под отдельные задачи, но не более того.
А вот структура catalog получилась значительно сложнее. Выглядит это вот так:
- catalog/:
- - io.py
- - funcs.py
- - catalog.py
- - pulsar.py
- - models.py

Где в io лежат по сути две функции load_data и dump_data для работы с pkl файлами + ещё несколько дополнительных функций чтобы доставать определенные файлы из папки .storage - например каталог ATNF как pandas DataFrame. Выглядит это примерно так:
```python
def load_data(path, name):
    data = pd.read_pickle(path / f"{name}.pkl")
    return data

def dump_data(data, path, name='data'):
    data.to_pickle(path / f"{name}.pkl")

# ATNF working

def load_atnf(path):
    storage = path.parent / ".storage"
    atnf = load_data(storage, 'atnf')
    return atnf

def download_atnf(path):
    atnf = psrqpy.QueryATNF().pandas
    dump_data(atnf, path, 'atnf')
    return atnf
```
В файле funcs.py у меня лежат три типа функций - переводы систем координат, создание массивов для разных распределений, и создание массива имен по кооординатам (что тоже можно отнести к системам координат). Это всё вспомогательные файлы, три главных по сути это catalog.py, pulsar.py и models.py. В pulsar.py у меня лежит вот такой код:
```python
import pandas as pd
from typing import NamedTuple

class Pulsar(NamedTuple):
    name: str
    ra: float # deg
    dec: float # deg
    f0: float = 1.0 # Hz
    f1: float = float('nan') # s^-2
    pmra: float = float('nan') # mas yr^-1
    pmdec: float = float('nan') # mas yr^-1
    px: float = float('nan') # mas
    dm: float = float('nan') # cm^-3 pc

atnf_format = {
    'name': 'PSRJ',
    'ra': 'RAJD',
    'dec': 'DECJD',
    'f0': 'F0',
    'f1': 'F1',
    'pmra': 'PMRA',
    'pmdec': 'PMDEC',
    'px': 'PX',
    'dm': 'DM',
}

def make_pulsar(params):
    return Pulsar(**params)

def make_pulsars(data):
    pulsars = []
    for row in data.itertuples():
        params = {attr: getattr(row, col) for attr, col in atnf_format.items()}
        pulsar = make_pulsar(params)
        pulsars.append(pulsar)
    return pulsars

def make_data(pulsars):
    data = pd.DataFrame(pulsars).rename(columns=atnf_format)
    return data
```
То есть, это просто абстракция на пульсар как набор параметров, словарь соответствия этих параметров названиям в каталоге ATNF и функции создания набора пульсаров по DataFrame и DataFrame по набору пульсаров. 

В файле catalog.py у меня лежит вот такой код:
```python
tempo_format = {
    'RAJD': 'RAJ',
    'DECJD': 'DECJ',
}

class Catalog:
    def __init__(self, data=None, name='test'):
        self.name = name
        required_cols = list(atnf_format.values()) 
        
        if data is None:
            data = pd.DataFrame(columns=required_cols)
        else:
            missing_cols = [col for col in required_cols if col not in data.columns]
            for col in missing_cols:
                if col == 'PSRJ':
                    data[col] = 'S' + data.index.astype(str)
                else:
                    data[col] = float('nan')
            
        self.data = data[required_cols]

    def pulsars(self):
        return make_pulsars(self.data)

    def tempo(self):
        return self.data.rename(columns=tempo_format)

    def add(self, pulsars):
        self.data = pd.concat([self.data, make_data(pulsars)], ignore_index=True)

    def sample(self, n, seed=42):
        sample_data = self.data.sample(n, random_state=seed)
        name =  f"{self.name}-n{n}-s{seed}"
        return Catalog(data=sample_data, name=name)

    def __repr__(self):
        return f"{self.name} catalog\n" + repr(self.data)

    def __str__(self):
        return f"{self.name} catalog\n" + str(self.data)

def make_catalog(n_psr, name='sample', method='test', params={}):
    data = make_synthetics(n_psr, method, params)
    return Catalog(data, name=name)

def load_catalog(path, name, prefix=True):
    filename = name
    if prefix:
        filename += "_cat"
    data = load_data(path, name=filename)
    return Catalog(data, name)

def save_catalog(catalog, path, prefix=True):
    filename = catalog.name
    if prefix:
        filename += "_cat"
    dump_data(catalog.data, path, name=filename)
```

То есть, это небольшая абстракция на DataFrame, чтобы с ней было просто удобнее работать - создавать DataFrame в нужном формате (с соответствующими Pulsar колонками), делать из него набор пульсаров и три функции для работы с Catalog - load, save и make. Load и Save просто сводят DataFrame к формату Catalog, а make_catalog использует модели синтетических данных для создания каталога. И тут мы переходим к третьему файлу - models.py, где по сути лежат много функций для генерации данных, функция   make_synthetics и словарь доступных методов. Выглядит это примерно так:
```python
def simple_ring(n_psr, seed_psr=42, alpha=10.0, ra_0=0.0, dec_0=0.0, radius=np.nan):
    phi_0, theta_0, _ = get_sph_coord(ra_0, dec_0)
    Phi, Theta, Rho = isotropic_ring(n_psr, seed_psr, radius, phi_0, theta_0, np.deg2rad(alpha))
    Ra, Dec, Px = sph2psr(Phi, Theta, Rho)
    Name = get_names(Ra, Dec, prefix='S')
    return [
        {'name': name, 'ra': ra, 'dec': dec, 'px': px} 
        for name, ra, dec, px in zip(Name, Ra, Dec, Px)
    ]

def make_synthetics(n_psr, method, params):
    pulsars = [] 
    synth_data = methods[method](n_psr, **params)
    for k in range(n_psr):
        psr = Pulsar(**synth_data[k])
        pulsars.append(psr)
    data = make_data(pulsars)
    return data 

methods = {
    'test': simple_test,
    'sphere': simple_sphere,
    'ball': simple_ball,
    'cap': simple_cap,
    'ring': simple_ring,
    'cone': simple_cone,
}
``` 

Вот, это на данный момент то, как выглядит мой код. У меня есть несколько фундаментальных вопросов по нему, которые я попробую аккуратно сформулировать. 
1) В идеале, я хочу добиться того, чтобы каждый модуль библиотеки ptatoolbox был независимым. То есть, у меня есть модули для разных задач, и они могут взаимодействовать, но только в рамках заданного API. То есть, условно говоря, модуль visuals умеет отрисовывать класс catalog, но он знает о нем только то, что к catalog можно обратиться по атрибуту, не больше. А Catalog работает со ссылками от DataManager, но ему нужна только ссылка в определенном формате, не больше. 

Поэтому мой первый вопрос звучит так - насколько хорошо прямо сейчас выглядит архитектура, взаимодействие частей кода между собой и т.д.? Может быть ты можешь предложить, как это ещё более аккуратно переписать?

2) Следующий шаг, который я хочу добавить в модуль catalog - это создание каталога на основании реальных данных с примесью искусственных. Условно говоря, с одной стороны у меня есть возможность сделать Catalog на основании DataFrame ATNF. С другой, я могу сделать синтетический Catalog с помощью make_synthetic, например пока что с искусственным распределением координат пульсаров по небу. А мне для моей задачи нужен формат, когда координаты генерируются искусственно, а остальные параметры случайно берутся из ATNF. Я понимаю, что я могу уже сейчас сгененрировать синтетический каталог с нужными координатами, и просто добавить в него данные из ATNF на уровне DataFrame, но я хочу это как-то обобщить до операций на уровне каталогов. Например, пока у меня есть такая мысль - добавить метод mix(self, real_catalog, paramas), который будет добавлять в свои данные, данные из real_catalog, и на его основании создавать новый каталог. Однако, как мне кажется это немного ломает идею абстрактной функции make_synthetic, в которую я просто задаю тип каталога и получаю нужный. Но с другой стороны, если делать это все на уровне make_synthetic, тогда нужно в качестве параметра опять таки давать реальные данные и значит опять у нас make_synthetic становится зависимый от импортов данных. 

Но есть например третий вариант - что мы сначала создаем реальный каталог на основании данных ATNF, а потом уже в него mix из данных синтетического каталога. Или как-то сделать этот препроцессинг при инициализации каталога. Не знаю. В общем, хотелось бы услышать твое мнение, как это правильно стоит делать.
3) Собственно, третья мысль которую хочется зафиксировать. Я хочу, чтобы мой класс Catalog стал универсальным форматом, необходимым и достаточным для анализа. То есть, мой класс Catalog - это корректная репрезентация массива пульсаров, и я работаю именно с каталогами, а не просто с данными. И вот хочется, чтобы на уровне реализации кода осталась эта идея, что я в конечном итоге пишу абстракцию Catalog, которая необходима и достаточна для работы с набором пульсаров с различными параметрами. 

Может быть в этом ключе ты мне тоже что-то подскажешь по архитектуре, как это можно улучшить.


