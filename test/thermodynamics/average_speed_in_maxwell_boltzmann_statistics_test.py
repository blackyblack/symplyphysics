from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import average_speed_in_maxwell_boltzmann_statistics

# Description
## The average speed of Argon (m = 39.948 u) at equilibrium temperature T = 100 K amounts
## to <v> = 230 m/s.

Args = namedtuple("Args", "t m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(100 * units.kelvin)
    m = Quantity(39.948 * units.amu)
    return Args(t=t, m=m)


def test_law(test_args: Args) -> None:
    result = average_speed_in_maxwell_boltzmann_statistics.calculate_average_speed(
        test_args.t, test_args.m)
    assert_equal(result, 230 * units.meter / units.second)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        average_speed_in_maxwell_boltzmann_statistics.calculate_average_speed(tb, test_args.m)
    with raises(TypeError):
        average_speed_in_maxwell_boltzmann_statistics.calculate_average_speed(100, test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        average_speed_in_maxwell_boltzmann_statistics.calculate_average_speed(test_args.t, mb)
    with raises(TypeError):
        average_speed_in_maxwell_boltzmann_statistics.calculate_average_speed(test_args.t, 100)
