import json

from aqbt.aquarium import serializer


def test_sample_serializer(aquarium):
    samples = aquarium.Sample.last(10)
    for sample in samples:
        data = serializer.sample_serializer(sample)
        print(json.dumps(data, indent=2))


def test_samples_serializer(aquarium):
    samples = aquarium.Sample.last(10)
    data = []
    with aquarium.with_cache(True) as sess:
        g = sess.browser.sample_network(samples)
        for model_name, sample_id in g.nodes:
            if model_name == "Sample":
                data.append(serializer.sample_serializer(sess.Sample.find(sample_id)))
    print(json.dumps(data, indent=2))
