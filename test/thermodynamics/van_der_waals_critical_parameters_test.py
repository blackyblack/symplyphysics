from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import van_der_waals_critical_parameters as parameters_law

# Description
## The critical parameters of Argon within the van der Waals equation of state model are
## V_c = 96 cm**3, p_c = 48.3 atm, T_c = 151 K. Van der Waals equation parameters for argon
## are a = 1.355 bar*(L/mol)**2 and b = 0.03201 L/mol.

Args = namedtuple("Args", "a b")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(1.355 * units.bar * (units.liter / units.mole)**2)
    b = Quantity(0.03201 * units.liter / units.mole)
    return Args(a=a, b=b)


def test_law(test_args: Args) -> None:
    Vc, pc, Tc = parameters_law.calculate_critical_parameters(test_args.a, test_args.b)
    assert_equal(Vc, 96 * units.centimeter**3 / units.mole)
    assert_equal(pc, 48.3 * units.atmosphere)
    assert_equal(Tc, 151 * units.kelvin)


def test_bad_first_parameter(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        parameters_law.calculate_critical_parameters(ab, test_args.b)
    with raises(TypeError):
        parameters_law.calculate_critical_parameters(100, test_args.b)


def test_bad_second_parameter(test_args: Args) -> None:
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        parameters_law.calculate_critical_parameters(test_args.a, bb)
    with raises(TypeError):
        parameters_law.calculate_critical_parameters(test_args.a, 100)
