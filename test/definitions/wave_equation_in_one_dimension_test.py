from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import wave_equation_in_one_dimension as wave_eqn

# Description
## A wave is traveling in a stretched string with a speed of 10 m/s. The position
## of the points of the string along the y-axis relative to the string's initial
## position is the displacement function in question. Assuming the solution of the
## wave equation takes the form of a cosine, the y-coordinate of a point of the wave
## at x = 0.3 m and t = 0.5 s is u = 5.54 cm, assuming the amplitude of the wave
## u_m = 10 cm and the phase lag is zero.

Args = namedtuple("Args", "u_m v phi x t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    u_m = Quantity(10.0 * units.centimeter)
    v = Quantity(10.0 * units.meter / units.second)
    phi = 0.0
    x = Quantity(0.3 * units.meter)
    t = Quantity(0.5 * units.second)
    return Args(u_m=u_m, v=v, phi=phi, x=x, t=t)


def test_definition(test_args: Args) -> None:
    result = wave_eqn.calculate_displacement(test_args.u_m, test_args.v, test_args.phi, test_args.x, test_args.t)
    assert_equal(result, 5.54 * units.centimeter)


def test_bad_phase_speed(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wave_eqn.calculate_displacement(test_args.u_m, vb, test_args.phi, test_args.x, test_args.t)
    with raises(TypeError):
        wave_eqn.calculate_displacement(test_args.u_m, 100, test_args.phi, test_args.x, test_args.t)


def test_bad_position(test_args: Args) -> None:
    xb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wave_eqn.calculate_displacement(test_args.u_m, test_args.v, test_args.phi, xb, test_args.t)
    with raises(TypeError):
        wave_eqn.calculate_displacement(test_args.u_m, test_args.v, test_args.phi, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        wave_eqn.calculate_displacement(test_args.u_m, test_args.v, test_args.phi, test_args.x, tb)
    with raises(TypeError):
        wave_eqn.calculate_displacement(test_args.u_m, test_args.v, test_args.phi, test_args.x, 100)
