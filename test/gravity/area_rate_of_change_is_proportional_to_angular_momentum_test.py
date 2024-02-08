from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity import (
    area_rate_of_change_is_proportional_to_angular_momentum as keplers_second_law,
)

# Description
## Earth's mass is 6.0e24 kg and its angular momentum is 2.7e40 kg*m**2/s.
## The rate of change of the area swept by it is about 2.25e9 km**2/s.


@fixture(name="test_args")
def test_args_fixture():
    L = Quantity(2.7e40 * units.kilogram * units.meter**2 / units.second)
    m = Quantity(6.0e24 * units.kilogram)
    Args = namedtuple("Args", "L m")
    return Args(L=L, m=m)


def test_law(test_args):
    result = keplers_second_law.calculate_rate_of_change_of_area(test_args.L, test_args.m)
    assert_equal(result, 2.25e9 * units.kilometer**2 / units.second)


def test_bad_angular_momentum(test_args):
    Lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        keplers_second_law.calculate_rate_of_change_of_area(Lb, test_args.m)
    with raises(TypeError):
        keplers_second_law.calculate_rate_of_change_of_area(100, test_args.m)


def test_bad_mass(test_args):
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        keplers_second_law.calculate_rate_of_change_of_area(test_args.L, mb)
    with raises(TypeError):
        keplers_second_law.calculate_rate_of_change_of_area(test_args.L, 100)
