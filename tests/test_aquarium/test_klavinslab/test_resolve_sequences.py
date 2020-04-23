from aqbt.aquarium.resolve_sequences import Resolver


def test_resolve_plasmids():
    session = sessions.aqproduction.with_cache()

    resolver = Resolver(sessions.aqproduction, sessions.klregistry)

    plasmids = session.Sample.last(
        10, query={"sample_type_id": resolver.plasmid_type.id}
    )

    for p in plasmids:
        resolver.resolve_sequence(p)


def test_resolve_plasmid():
    resolver = Resolver(sessions.aqproduction, sessions.klregistry)
    resolver.force_build_at_depths = [0, 1]
    resolver.resolve_sequence(sessions.aqproduction.Sample.find(32869))


def test_resolve_fragment():
    resolver = Resolver(sessions.aqproduction, sessions.klregistry)
    resolver.force_build_at_depths = [0]
    seq = resolver.resolve_sequence(sessions.aqproduction.Sample.find(32868))

    for a in seq.annotations:
        print(a.name, a.start, a.end, len(seq.bases))


def test_unregistered():
    resolver = Resolver(sessions.aqproduction, sessions.klregistry)
    samples = sessions.aqproduction.Sample.last(
        10, query={"sample_type_id": [d.id for d in resolver.dna_types]}
    )
    for s in samples:
        resolver.force_build_at_depths = [0, 1]
        seq = resolver.resolve_sequence(s)
