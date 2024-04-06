from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import specific_heat_capacity

# Description
## The specific heat capacity of 2 kg of a substance with heat capacity C = 10 J/K is c = 5 J/(K*kg).

Args = namedtuple("Args", "c m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    c = Quantity(10 * units.joule / units.kelvin)
    m = Quantity(2 * units.kilogram)
    return Args(c=c, m=m)


def test_law(test_args: Args) -> None:
    result = specific_heat_capacity.calculate_specific_heat_capacity(test_args.c, test_args.m)
    assert_equal(result, 5 * units.joule / (units.kelvin * units.kilogram))


def test_bad_heat_capacity(test_args: Args) -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        specific_heat_capacity.calculate_specific_heat_capacity(cb, test_args.m)
    with raises(TypeError):
        specific_heat_capacity.calculate_specific_heat_capacity(100, test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        specific_heat_capacity.calculate_specific_heat_capacity(test_args.c, mb)
    with raises(TypeError):
        specific_heat_capacity.calculate_specific_heat_capacity(test_args.c, 100)
