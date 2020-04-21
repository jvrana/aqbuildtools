from aqbt import utils
from colour import Color


def test_random_color():
    for _ in range(1000):
        c = utils.random_color()
        print(c)
        assert len(c) == 7
        Color(c)
