from collections import namedtuple
from pytest import fixture
from symplyphysics import (
    assert_equal,
    units,
)
from symplyphysics.laws.waves import wavespeed_from_medium as speed_law

# Description
## Known propagation speed of electromagnetical wave in air (with it's refraction factor of 1.0003) is 2.997925e8 km/s.


@fixture(name="test_args")
def test_args_fixture():
    refraction_factor = 1.0003
    Args = namedtuple("Args", ["refraction_factor"])
    return Args(refraction_factor=refraction_factor)


def test_basic_speed(test_args):
    result = speed_law.calculate_wavespeed(test_args.refraction_factor)
    assert_equal(result, 299792500 * units.meter / units.second)
