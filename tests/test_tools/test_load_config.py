from aqbt import AquariumBuildTools
from aqbt.tools import parse_config
import toml


def test_parse_config(config_path):
    config = toml.load(config_path)
    config = parse_config(config)
    assert config['default']['aquarium']['url'] == 'http://0.0.0.0/'
    assert config['local']['aquarium']['url'] == 'http://10.10.0.0/'


def test_aqbt_init(config_path):
    aqbt = AquariumBuildTools.from_toml(config_path)
