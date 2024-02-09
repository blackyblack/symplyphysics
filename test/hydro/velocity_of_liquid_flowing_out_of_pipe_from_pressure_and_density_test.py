from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.hydro import velocity_of_liquid_flowing_out_of_pipe_from_pressure_and_density as velocity_law

# Description
## The density of water is 997 [kilogram / meter^3]. At a pressure of 5e5 pascal,
## the water velocity will be 31.67 meter per second.
## https://www.indigomath.ru//raschety/C7zPVP.html

Args = namedtuple("Args", ["pressure", "density"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    pressure = Quantity(5e5 * units.pascal)
    density = Quantity(997 * (units.kilogram / units.meter**3))
    return Args(pressure=pressure, density=density)


def test_basic_velocity(test_args: Args) -> None:
    result = velocity_law.calculate_velocity(test_args.pressure, test_args.density)
    assert_equal(result, 31.67 * units.meter / units.second)


def test_bad_pressure(test_args: Args) -> None:
    pressure = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        velocity_law.calculate_velocity(pressure, test_args.density)
    with raises(TypeError):
        velocity_law.calculate_velocity(100, test_args.density)


def test_bad_density(test_args: Args) -> None:
    density = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        velocity_law.calculate_velocity(test_args.pressure, density)
    with raises(TypeError):
        velocity_law.calculate_velocity(test_args.pressure, 100)
