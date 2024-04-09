from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import critical_pressure as critical_law

# Description
## The critical pressure of Argon within the van der Waals equation of state model is
## p_c = 48.3 atm. Van der Waals equation parameters for argon are a = 1.355 bar*(L/mol)**2
## and b = 0.03201 L/mol.

Args = namedtuple("Args", "a b")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(1.355 * units.bar * (units.liter / units.mole)**2)
    b = Quantity(0.03201 * units.liter / units.mole)
    return Args(a=a, b=b)


def test_law(test_args: Args) -> None:
    pc = critical_law.calculate_critical_pressure(test_args.a, test_args.b)
    assert_equal(pc, 48.3 * units.atmosphere)


def test_bad_first_parameter(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        critical_law.calculate_critical_pressure(ab, test_args.b)


def test_bad_second_parameter(test_args: Args) -> None:
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        critical_law.calculate_critical_pressure(test_args.a, bb)
