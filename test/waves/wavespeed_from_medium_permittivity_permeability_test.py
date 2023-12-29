from collections import namedtuple
from pytest import approx, fixture
from symplyphysics import (
    units,
    SI,
    convert_to,
)
from symplyphysics.laws.waves import wavespeed_from_medium_permittivity_permeability as speed_law

# Description
## Known propagation speed of electromagnetical wave in vacuum is 299792458 m/s.


@fixture(name="test_args")
def test_args_fixture():
    relative_permittivity = 1
    relative_permeability = 1
    Args = namedtuple("Args", ["relative_permittivity", "relative_permeability"])
    return Args(relative_permittivity=relative_permittivity,
                relative_permeability=relative_permeability)


def test_basic_speed(test_args):
    result = speed_law.calculate_wavespeed(
        test_args.relative_permittivity,
        test_args.relative_permeability
    )
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.speed)
    result_freq_1 = convert_to(result, units.meter / units.second).evalf(5)
    assert result_freq_1 == approx(299792458, rel=1e-2)
