from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics.equations_of_state import dieterici_equation

# Description
## For a fluid described by the Dieterici equation of molar volume V = 0.3 m**3/mol,
## temperature T = 200 K, and parameters a = 2 bar*(L/mol)**2 and b = 0.1 L/mol,
## the pressure is 5.54 kPa.

Args = namedtuple("Args", "vm t a b")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    vm = Quantity(0.3 * units.meter**3 / units.mole)
    t = Quantity(200 * units.kelvin)
    a = Quantity(2 * units.bar * (units.liter / units.mole)**2)
    b = Quantity(0.1 * units.liter / units.mole)
    return Args(vm=vm, t=t, a=a, b=b)


def test_law(test_args: Args) -> None:
    result = dieterici_equation.calculate_pressure(test_args.vm, test_args.t, test_args.a,
        test_args.b)
    assert_equal(result, 5.54 * prefixes.kilo * units.pascal)


def test_bad_molar_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        dieterici_equation.calculate_pressure(vb, test_args.t, test_args.a, test_args.b)
    with raises(TypeError):
        dieterici_equation.calculate_pressure(100, test_args.t, test_args.a, test_args.b)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        dieterici_equation.calculate_pressure(test_args.vm, tb, test_args.a, test_args.b)
    with raises(TypeError):
        dieterici_equation.calculate_pressure(test_args.vm, 100, test_args.a, test_args.b)


def test_bad_first_parameter(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        dieterici_equation.calculate_pressure(test_args.vm, test_args.t, ab, test_args.b)
    with raises(TypeError):
        dieterici_equation.calculate_pressure(test_args.vm, test_args.t, 100, test_args.b)


def test_bad_second_parameter(test_args: Args) -> None:
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        dieterici_equation.calculate_pressure(test_args.vm, test_args.t, test_args.a, bb)
    with raises(TypeError):
        dieterici_equation.calculate_pressure(test_args.vm, test_args.t, test_args.a, 100)
