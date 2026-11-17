import ptatoolbox as pta
import typer

def main(exp: str = 'test-io'):
    dm = pta.DataManager()
    print(dm.root)
    print(dm.data)
    print(dm.config)
    print(dm.raw)
    print(dm.processed)

    sim_path = dm.make_dir(path=dm.raw, name=exp)
    res_path = dm.make_dir(path=dm.processed, name=exp)
    stat_path = dm.make_dir(path=res_path, name='stat')

    cfg_path = dm.config / f"{exp}.yaml"

    print(sim_path)
    print(res_path)
    print(stat_path)
    print(cfg_path)
    
    config = pta.parse_yaml(path=cfg_path)
    print(config)

    return None

if __name__ == '__main__':
    typer.run(main)