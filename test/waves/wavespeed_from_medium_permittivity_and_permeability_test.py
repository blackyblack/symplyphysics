from collections import namedtuple
from pytest import fixture
from symplyphysics import (
    assert_equal,
    units,
)
from symplyphysics.laws.waves import wavespeed_from_medium_permittivity_and_permeability as speed_law

# Description
## Known propagation speed of electromagnetical wave in vacuum is 299792458 m/s.

Args = namedtuple("Args", ["relative_permittivity", "relative_permeability"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 1
    relative_permeability = 1
    return Args(relative_permittivity=relative_permittivity,
        relative_permeability=relative_permeability)


def test_basic_speed(test_args: Args) -> None:
    result = speed_law.calculate_wavespeed(test_args.relative_permittivity,
        test_args.relative_permeability)
    assert_equal(result, 299792458 * units.meter / units.second)
