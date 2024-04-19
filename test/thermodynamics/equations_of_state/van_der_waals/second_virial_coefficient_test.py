from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import (
    second_virial_coefficient as virial_law,)

# Description
## The van der Waals coefficients for Argon are a = 1.355 bar*(L/mol)**2 and
## b = 0.032 L/mol. At temperature T = 300 K the second virial coefficient in the
## model of van der Waals equation of state amounts to -22.3 mL/mol.

Args = namedtuple("Args", "a b t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(1.355 * units.bar * (units.liter / units.mole)**2)
    b = Quantity(0.032 * units.liter / units.mole)
    t = Quantity(300 * units.kelvin)
    return Args(a=a, b=b, t=t)


def test_law(test_args: Args) -> None:
    result = virial_law.calculate_second_virial_coefficient(test_args.a, test_args.b, test_args.t)
    assert_equal(result, -22.3 * units.milliliter / units.mole, tolerance=2e-3)


def test_bad_first_coefficient(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        virial_law.calculate_second_virial_coefficient(ab, test_args.b, test_args.t)
    with raises(TypeError):
        virial_law.calculate_second_virial_coefficient(100, test_args.b, test_args.t)


def test_bad_second_coefficient(test_args: Args) -> None:
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        virial_law.calculate_second_virial_coefficient(test_args.a, bb, test_args.t)
    with raises(TypeError):
        virial_law.calculate_second_virial_coefficient(test_args.a, 100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        virial_law.calculate_second_virial_coefficient(test_args.a, test_args.b, tb)
    with raises(TypeError):
        virial_law.calculate_second_virial_coefficient(test_args.a, test_args.b, 100)
