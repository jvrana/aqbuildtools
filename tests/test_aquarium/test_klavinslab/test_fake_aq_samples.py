from aqbt.aquarium.faker import FakeSampleGenerator

def test_fake_primer(aquarium):
    gen = FakeSampleGenerator(aquarium)
    sample = gen.fake_primer()
    print(sample.properties)


def test_fake_plasmid(aquarium):
    gen = FakeSampleGenerator(aquarium)
    sample = gen.fake_plasmid()
    print(sample.properties)


def test_fake_fragment(aquarium):
    gen = FakeSampleGenerator(aquarium)
    sample = gen.fake_fragment()


def test_fake_library(aquarium):
    gen = FakeSampleGenerator(aquarium)
    assert gen.make_fake_library(50, 50, 50)
