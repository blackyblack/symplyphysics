from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import heat_capacity_via_molar_heat_capacity

# Description
## The heat capacity of 2 mol of a substance with molar heat capacity C_m = 5 J/(K*mol) is C = 10 J/K.

Args = namedtuple("Args", "c n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    c = Quantity(5 * units.joule / (units.kelvin * units.mole))
    n = Quantity(2 * units.mole)
    return Args(c=c, n=n)


def test_law(test_args: Args) -> None:
    result = heat_capacity_via_molar_heat_capacity.calculate_heat_capacity(test_args.c, test_args.n)
    assert_equal(result, 10 * units.joule / units.kelvin)


def test_bad_molar_heat_capacity(test_args: Args) -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heat_capacity_via_molar_heat_capacity.calculate_heat_capacity(cb, test_args.n)
    with raises(TypeError):
        heat_capacity_via_molar_heat_capacity.calculate_heat_capacity(100, test_args.n)


def test_bad_amount_of_substance(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heat_capacity_via_molar_heat_capacity.calculate_heat_capacity(test_args.c, nb)
    with raises(TypeError):
        heat_capacity_via_molar_heat_capacity.calculate_heat_capacity(test_args.c, 100)
