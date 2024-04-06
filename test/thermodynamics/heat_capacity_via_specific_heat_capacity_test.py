from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import heat_capacity_via_specific_heat_capacity

# Description
## The heat capacity of 2 kg of a substance with specific heat capacity c = 5 J/(K*kg) is C = 10 J/K.

Args = namedtuple("Args", "c m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    c = Quantity(5 * units.joule / (units.kelvin * units.kilogram))
    m = Quantity(2 * units.kilogram)
    return Args(c=c, m=m)


def test_law(test_args: Args) -> None:
    result = heat_capacity_via_specific_heat_capacity.calculate_heat_capacity(test_args.c, test_args.m)
    assert_equal(result, 10 * units.joule / units.kelvin)


def test_bad_specific_heat_capacity(test_args: Args) -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heat_capacity_via_specific_heat_capacity.calculate_heat_capacity(cb, test_args.m)
    with raises(TypeError):
        heat_capacity_via_specific_heat_capacity.calculate_heat_capacity(100, test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heat_capacity_via_specific_heat_capacity.calculate_heat_capacity(test_args.c, mb)
    with raises(TypeError):
        heat_capacity_via_specific_heat_capacity.calculate_heat_capacity(test_args.c, 100)
