from aqbt.aquarium.resolve_sequences import Resolver


def test_resolve_plasmids(aquarium, registry):
    session = aquarium.with_cache()

    resolver = Resolver(session, registry)

    plasmids = session.Sample.last(
        10, query={"sample_type_id": resolver.plasmid_type.id}
    )

    for p in plasmids:
        resolver.resolve_sequence(p)


def test_resolve_plasmid(aquarium, registry):
    resolver = Resolver(aquarium, registry)
    resolver.force_build_at_depths = [0, 1]
    resolver.resolve_sequence(aquarium.Sample.find(32869))


def test_resolve_fragment(aquarium, registry):
    resolver = Resolver(aquarium, registry)
    resolver.force_build_at_depths = [0]
    seq = resolver.resolve_sequence(aquarium.Sample.find(32868))

    for a in seq.annotations:
        print(a.name, a.start, a.end, len(seq.bases))


def test_unregistered(aquarium, registry):
    resolver = Resolver(aquarium, registry)
    samples = aquarium.Sample.last(
        10, query={"sample_type_id": [d.id for d in resolver.dna_types]}
    )
    for s in samples:
        resolver.force_build_at_depths = [0, 1]
        seq = resolver.resolve_sequence(s)
