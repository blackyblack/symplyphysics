from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import braking_path as path_law

# Description
## Let the mass of the object be 10 kilogram, the speed is 20 meter per second, and the friction force 
## is 40 newton. Then the stopping distance will be 50 meter.
## https://www.indigomath.ru/formuly-po-fizike/dinamika/put-tormozhenija-2.html


@fixture(name="test_args")
def test_args_fixture():
    mass = Quantity(10 * units.kilogram)
    velocity = Quantity(20 * units.meter / units.second)
    friction_force = Quantity(40 * units.newton)

    Args = namedtuple("Args", ["mass", "velocity", "friction_force"])
    return Args(mass=mass, velocity=velocity, friction_force=friction_force)


def test_basic_braking_path(test_args):
    result = path_law.calculate_braking_path(test_args.mass, test_args.velocity, test_args.friction_force)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_voltage = convert_to(result, units.meter).evalf(5)
    assert result_voltage == approx(50, 0.01)


def test_bad_mass(test_args):
    mass = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        path_law.calculate_braking_path(mass, test_args.velocity, test_args.friction_force)
    with raises(TypeError):
        path_law.calculate_braking_path(100, test_args.velocity, test_args.friction_force)


def test_bad_velocity(test_args):
    velocity = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        path_law.calculate_braking_path(test_args.mass, velocity, test_args.friction_force)
    with raises(TypeError):
        path_law.calculate_braking_path(test_args.mass, 100, test_args.friction_force)


def test_bad_friction_force(test_args):
    friction_force = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        path_law.calculate_braking_path(test_args.mass, test_args.velocity, friction_force)
    with raises(TypeError):
        path_law.calculate_braking_path(test_args.mass, test_args.velocity, 100)
