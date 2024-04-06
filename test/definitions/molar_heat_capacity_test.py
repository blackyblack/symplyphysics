from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import molar_heat_capacity

# Description
## The molar heat capacity of 2 mol of a substance with heat capacity C = 10 J/K is C_m = 5 J/(k*mol).

Args = namedtuple("Args", "c n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    c = Quantity(10 * units.joule / units.kelvin)
    n = Quantity(2 * units.mole)
    return Args(c=c, n=n)


def test_law(test_args: Args) -> None:
    result = molar_heat_capacity.calculate_molar_heat_capacity(test_args.c, test_args.n)
    assert_equal(result, 5 * units.joule / (units.kelvin * units.mole))


def test_bad_heat_capacity(test_args: Args) -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        molar_heat_capacity.calculate_molar_heat_capacity(cb, test_args.n)
    with raises(TypeError):
        molar_heat_capacity.calculate_molar_heat_capacity(100, test_args.n)


def test_bad_amount_of_substance(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        molar_heat_capacity.calculate_molar_heat_capacity(test_args.c, nb)
    with raises(TypeError):
        molar_heat_capacity.calculate_molar_heat_capacity(test_args.c, 100)
