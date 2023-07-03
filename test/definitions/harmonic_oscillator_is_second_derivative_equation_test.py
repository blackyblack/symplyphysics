from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import pi
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import harmonic_oscillator_is_second_derivative_equation as oscillator

# Description
## Oscillator has an angular frequency of 2*pi herts. Therefore it's linear frequency is 1 herts. After
## 1.5 seconds it will be on the negative side of the equilibrium, i.e. at -20 meters displacement.


@fixture(name="test_args")
def test_args_fixture():
    A = Quantity(20 * units.meter)
    w = Quantity(2 * pi / 1 * units.hertz)
    t = Quantity(1.5 * units.second)
    Args = namedtuple("Args", ["A", "w", "t"])
    return Args(A=A, w=w, t=t)


def test_basic_displacement(test_args):
    result = oscillator.calculate_displacement(test_args.A, test_args.w, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_frequency = convert_to(result, units.meter).evalf(4)
    assert result_frequency == approx(-20, 0.01)


def test_bad_amplitude(test_args):
    Ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        oscillator.calculate_displacement(Ab, test_args.w, test_args.t)
    with raises(TypeError):
        oscillator.calculate_displacement(100, test_args.w, test_args.t)


def test_bad_frequency(test_args):
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        oscillator.calculate_displacement(test_args.A, wb, test_args.t)
    with raises(TypeError):
        oscillator.calculate_displacement(test_args.A, 100, test_args.t)


def test_bad_time(test_args):
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        oscillator.calculate_displacement(test_args.A, test_args.w, tb)
    with raises(TypeError):
        oscillator.calculate_displacement(test_args.A, test_args.w, 100)
