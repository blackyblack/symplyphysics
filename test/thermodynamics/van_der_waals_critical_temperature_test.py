from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import van_der_waals_critical_temperature as critical_law

# Description
## The critical temperature of Argon within the van der Waals equation of state model are
## T_c = 151 K. Van der Waals equation parameters for argon are a = 1.355 bar*(L/mol)**2
## and b = 0.03201 L/mol.

Args = namedtuple("Args", "a b")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(1.355 * units.bar * (units.liter / units.mole)**2)
    b = Quantity(0.03201 * units.liter / units.mole)
    return Args(a=a, b=b)


def test_law(test_args: Args) -> None:
    tc = critical_law.calculate_critical_temperature(test_args.a, test_args.b)
    assert_equal(tc, 151 * units.kelvin)


def test_bad_first_parameter(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        critical_law.calculate_critical_temperature(ab, test_args.b)


def test_bad_second_parameter(test_args: Args) -> None:
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        critical_law.calculate_critical_temperature(test_args.a, bb)
