from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import braking_path as path_law

# Description
## Let the mass of the object be 10 kilogram, the speed is 20 meter per second, and the friction force
## is 40 newton. Then the stopping distance will be 50 meter.
## https://www.indigomath.ru/formuly-po-fizike/dinamika/put-tormozhenija-2.html

Args = namedtuple("Args", ["mass", "velocity", "friction_force"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mass = Quantity(10 * units.kilogram)
    velocity = Quantity(20 * units.meter / units.second)
    friction_force = Quantity(40 * units.newton)
    return Args(mass=mass, velocity=velocity, friction_force=friction_force)


def test_basic_braking_path(test_args: Args) -> None:
    result = path_law.calculate_braking_path(test_args.mass, test_args.velocity,
        test_args.friction_force)
    assert_equal(result, 50 * units.meter)


def test_bad_mass(test_args: Args) -> None:
    mass = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        path_law.calculate_braking_path(mass, test_args.velocity, test_args.friction_force)
    with raises(TypeError):
        path_law.calculate_braking_path(100, test_args.velocity, test_args.friction_force)


def test_bad_velocity(test_args: Args) -> None:
    velocity = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        path_law.calculate_braking_path(test_args.mass, velocity, test_args.friction_force)
    with raises(TypeError):
        path_law.calculate_braking_path(test_args.mass, 100, test_args.friction_force)


def test_bad_friction_force(test_args: Args) -> None:
    friction_force = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        path_law.calculate_braking_path(test_args.mass, test_args.velocity, friction_force)
    with raises(TypeError):
        path_law.calculate_braking_path(test_args.mass, test_args.velocity, 100)
