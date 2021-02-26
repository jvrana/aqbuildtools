import pytest

from aqbt.aquarium.registry import FakeRegistryConnector


def test_fake_registry_connector():

    registry = FakeRegistryConnector((1000, 5000), (100, 2000))
    seqs = registry.all(100)
    assert len(seqs) == 100
    print(seqs)
    print(seqs[0].keys())


def test_fake_registry_connector_find():

    registry = FakeRegistryConnector((1000, 5000), (100, 2000))
    seqs = registry.all(100)
    assert len(seqs) == 100

    id = seqs[0]["entityRegistryId"]
    print(id)
    seq = registry.find(id)
    assert seq


@pytest.mark.parametrize("always_return", [False, True])
def test_fake_registry_connector_find_missing(always_return):

    registry = FakeRegistryConnector((1000, 5000), (100, 2000))
    seqs = registry.all(100)
    assert len(seqs) == 100

    seq = registry.find("asdfasdf", always_return=always_return)
    assert (seq is None) != always_return
