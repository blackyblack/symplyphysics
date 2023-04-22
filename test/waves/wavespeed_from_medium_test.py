from collections import namedtuple
from pytest import approx, fixture
from symplyphysics import (
    units,
    SI,
    convert_to,
)
from symplyphysics.laws.waves import wavespeed_from_medium as speed_law

# Description
## Known propagation speed of electromagnetical wave in air (with it's refraction factor of 1.0003) is 2.997925e8 km/s.


@fixture
def test_args():
    refraction_factor = 1.0003
    Args = namedtuple("Args", ["refraction_factor"])
    return Args(refraction_factor=refraction_factor)


def test_basic_speed(test_args):
    result = speed_law.calculate_wavespeed(test_args.refraction_factor)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.speed)
    result_freq_1 = convert_to(result, units.meter / units.second).subs({
        units.meter: 1,
        units.second: 1
    }).evalf(5)
    assert result_freq_1 == approx(299792500, 1)
