from os.path import join, dirname, abspath

from BCBio import GFF
from aqbt.genome import GenomeIntegrator
import pytest
from aqbt import biopython



@pytest.fixture(scope="function")
def records(fixtures_path):
    return list(GFF.parse([join(fixtures_path, "cenpk.gff")]))


@pytest.mark.parametrize("locus", [(9500, 10000, 11000, 11500)])
@pytest.mark.parametrize("rc", [True, False], ids=["reverse_complement cassette", ""])
def test_basic(records, locus, rc):
    a, b, c, d = locus
    chromosomes = records[:]
    new_chr = biopython.random_record(20000)
    chromosomes.append(new_chr)
    cassette = biopython.random_record(5000, name="my_feature", auto_annotate=True)
    hom1 = new_chr[a:b]
    hom2 = new_chr[c:d]
    integration_cassette = hom1 + cassette + hom2
    if rc:
        integration_cassette = integration_cassette.reverse_complement()
    integrator = GenomeIntegrator()
    integrator.logger.set_level("DEBUG")
    biopython.make_linear([integration_cassette])
    biopython.make_linear(chromosomes)
    new_genome = integrator.integrate(chromosomes, [integration_cassette])

    expected_chr = new_chr[:b] + cassette + new_chr[c:]
    assert new_genome
    assert str(cassette.seq) in str(new_genome["genome"][-1].seq)
    assert str(cassette.seq) not in str(new_genome["genome"][-2].seq)
    assert str(new_chr[:b].seq) in str(new_genome["genome"][-1].seq)
    assert str(new_chr[c:].seq) in str(new_genome["genome"][-1].seq)
    assert str(new_genome["genome"][-1].seq) == str(expected_chr.seq)

    labels = [f.qualifiers["label"][0] for f in new_genome["genome"][-1].features]
    assert "my_feature" in labels


@pytest.mark.parametrize("locus", [(9500, 10000, 11000, 11500)])
@pytest.mark.parametrize("rc", [True, False], ids=["reverse_complement cassette", ""])
def test_basic_integrate_on_bottom_strand(records, locus, rc):
    a, b, c, d = locus
    chromosomes = records[:]
    new_chr = biopython.random_record(20000)
    chromosomes.append(new_chr)
    cassette = biopython.random_record(5000)
    hom1 = new_chr[a:b]
    hom2 = new_chr[c:d]
    integration_cassette = (
        hom2.reverse_complement() + cassette + hom1.reverse_complement()
    )
    if rc:
        integration_cassette = integration_cassette.reverse_complement()
    integrator = GenomeIntegrator()
    integrator.logger.set_level("DEBUG")
    biopython.make_linear([integration_cassette])
    biopython.make_linear(chromosomes)
    new_genome = integrator.integrate(chromosomes, [integration_cassette])

    expected_chr = new_chr[:b] + cassette.reverse_complement() + new_chr[c:]
    assert new_genome
    assert str(cassette.reverse_complement().seq) in str(new_genome["genome"][-1].seq)
    assert str(cassette.reverse_complement().seq) not in str(
        new_genome["genome"][-2].seq
    )
    assert str(new_chr[:b].seq) in str(new_genome["genome"][-1].seq)
    assert str(new_chr[c:].seq) in str(new_genome["genome"][-1].seq)
    assert str(new_genome["genome"][-1].seq) == str(expected_chr.seq)


@pytest.mark.parametrize("locus", [(11000, 11500, 9500, 10000)])
@pytest.mark.parametrize("rc", [True, False], ids=["reverse_complement cassette", ""])
def test_fail(records, locus, rc):
    a, b, c, d = locus
    chromosomes = records[:]
    new_chr = biopython.random_record(20000)
    chromosomes.append(new_chr)
    cassette = biopython.random_record(5000)
    hom1 = new_chr[a:b]
    hom2 = new_chr[c:d]
    integration_cassette = (
        hom2.reverse_complement() + cassette + hom1.reverse_complement()
    )
    if rc:
        integration_cassette = integration_cassette.reverse_complement()
    integrator = GenomeIntegrator()
    integrator.logger.set_level("DEBUG")
    biopython.make_linear([integration_cassette])
    biopython.make_linear(chromosomes)
    new_genome = integrator.integrate(chromosomes, [integration_cassette])
    assert new_genome["errors"]


@pytest.mark.parametrize("locus", [(11000, 11500, 9500, 10000)])
@pytest.mark.parametrize("rc", [True, False], ids=["reverse_complement cassette", ""])
def test_fail_integrate_on_bottom_strand(records, locus, rc):
    a, b, c, d = locus
    chromosomes = records[:]
    new_chr = biopython.random_record(20000)
    chromosomes.append(new_chr)
    cassette = biopython.random_record(5000)
    hom1 = new_chr[a:b]
    hom2 = new_chr[c:d]
    integration_cassette = (
        hom2.reverse_complement() + cassette + hom1.reverse_complement()
    )
    if rc:
        integration_cassette = integration_cassette.reverse_complement()
    integrator = GenomeIntegrator()
    integrator.logger.set_level("DEBUG")
    biopython.make_linear([integration_cassette])
    biopython.make_linear(chromosomes)
    new_genome = integrator.integrate(chromosomes, [integration_cassette])
    assert new_genome["errors"]