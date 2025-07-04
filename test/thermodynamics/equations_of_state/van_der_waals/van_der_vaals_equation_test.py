from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import van_der_vaals_equation as waals_law

# Description
## Test example from https://studfile.net/preview/1772224/page:8/ example 2.24

# Find pressure in state 1. For this example koefficients for van-der-vaals state equation

# a = (V/nu)^2 * (p2 * T1 - p1 * T2) / (T2 - T1) = 0.191 [Pa * m^6 / mole^2]
# b = V - nu * R * (T2 - T1) / (p2 - p1) = 4.532 * 10^(-5) [m^3 / mole]
# V = 0.25 [liters]
# nu = 1 [moles]
# T = 300 [K]

# Pressure should be 90 atm

Args = namedtuple("Args", ["t", "v", "a", "b"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(300 * units.kelvins)
    v_ = Quantity(0.25 * units.liters)
    nu = Quantity(1 * units.mole)
    v = Quantity(v_ / nu)
    a = Quantity(0.191 * units.pascals * (units.meter**3 / units.mole)**2)
    b = Quantity(4.532 * 1e-5 * units.meters**3 / units.moles)
    return Args(t=t, v=v, a=a, b=b)


def test_basic_pressure(test_args: Args) -> None:
    result = waals_law.calculate_pressure(test_args.v, test_args.t, test_args.a, test_args.b)
    assert_equal(result, 90 * units.atm, relative_tolerance=0.01)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        waals_law.calculate_pressure(test_args.v, tb, test_args.a, test_args.b)
    with raises(TypeError):
        waals_law.calculate_pressure(test_args.v, 100, test_args.a, test_args.b)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        waals_law.calculate_pressure(vb, test_args.t, test_args.a, test_args.b)
    with raises(TypeError):
        waals_law.calculate_pressure(100, test_args.t, test_args.a, test_args.b)


def test_bad_molecules_volume_parameter(test_args: Args) -> None:
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        waals_law.calculate_pressure(test_args.v, test_args.t, test_args.a, bb)
    with raises(TypeError):
        waals_law.calculate_pressure(test_args.v, test_args.t, test_args.a, 100)


def test_bad_bonding_forces_parameter(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        waals_law.calculate_pressure(test_args.v, test_args.t, ab, test_args.b)
    with raises(TypeError):
        waals_law.calculate_pressure(test_args.v, test_args.t, 100, test_args.b)
