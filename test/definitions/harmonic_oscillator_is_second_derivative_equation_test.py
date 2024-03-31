from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import harmonic_oscillator_is_second_derivative_equation as oscillator

# Description
## Oscillator has an angular frequency of 2*pi herts. Therefore it's linear frequency is 1 herts. After
## 1.5 seconds it will be on the negative side of the equilibrium, i.e. at -20 meters displacement.

Args = namedtuple("Args", ["A", "w", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    A = Quantity(20 * units.meter)
    w = Quantity(2 * pi * units.hertz)
    t = Quantity(1.5 * units.second)
    return Args(A=A, w=w, t=t)


def test_basic_displacement(test_args: Args) -> None:
    result = oscillator.calculate_displacement(test_args.A, test_args.w, test_args.t)
    assert_equal(result, -20 * units.meter)


def test_bad_amplitude(test_args: Args) -> None:
    Ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        oscillator.calculate_displacement(Ab, test_args.w, test_args.t)
    with raises(TypeError):
        oscillator.calculate_displacement(100, test_args.w, test_args.t)


def test_bad_frequency(test_args: Args) -> None:
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        oscillator.calculate_displacement(test_args.A, wb, test_args.t)
    with raises(TypeError):
        oscillator.calculate_displacement(test_args.A, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        oscillator.calculate_displacement(test_args.A, test_args.w, tb)
    with raises(TypeError):
        oscillator.calculate_displacement(test_args.A, test_args.w, 100)
