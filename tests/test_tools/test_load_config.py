from aqbt import AquariumBuildTools


def test_aqbt_init(config_path):
    aqbt = AquariumBuildTools.from_toml(config_path)