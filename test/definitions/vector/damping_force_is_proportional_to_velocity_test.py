from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, Quantity, errors
from symplyphysics.definitions.vector import (damping_force_is_proportional_to_velocity as
    damping_def)

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

# Description
## A damping force acts on a body. The body's velocity is (4, -2, 0) m/s, the damping constant
## of the system is 0.05 kg/s. Then the damping force is (-0.2, 0.1, 0.0) N.

Args = namedtuple("Args", "b v f")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    b = Quantity(0.05 * units.kilogram / units.second)
    v = QuantityCoordinateVector([
        Quantity(4.0 * units.meter / units.second),
        Quantity(-2.0 * units.meter / units.second),
        Quantity(0.0 * units.meter / units.second),
    ], CARTESIAN)
    f = QuantityCoordinateVector([
        Quantity(-0.2 * units.newton),
        Quantity(0.1 * units.newton),
        Quantity(0.0 * units.newton),
    ], CARTESIAN)
    return Args(b=b, v=v, f=f)


def test_damping_force_definition(test_args: Args) -> None:
    result = damping_def.calculate_damping_force(test_args.b, test_args.v)
    assert_equal_vectors(result, test_args.f)


def test_velocity_law(test_args: Args) -> None:
    result = damping_def.calculate_velocity(test_args.b, test_args.f)
    assert_equal_vectors(result, test_args.v)


def test_bad_damping_constant(test_args: Args) -> None:
    bb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        damping_def.calculate_damping_force(bb, test_args.v)
    with raises(TypeError):
        damping_def.calculate_damping_force(100, test_args.v)
    with raises(errors.UnitsError):
        damping_def.calculate_velocity(bb, test_args.f)
    with raises(TypeError):
        damping_def.calculate_velocity(100, test_args.f)


def test_bad_velocity(test_args: Args) -> None:
    v_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        damping_def.calculate_damping_force(test_args.b, v_bad_vector)

    v_scalar = Quantity(1.0 * units.meter / units.second)
    with raises(ValueError):
        damping_def.calculate_damping_force(test_args.b, v_scalar)

    with raises(TypeError):
        damping_def.calculate_damping_force(test_args.b, 100)
    with raises(TypeError):
        damping_def.calculate_damping_force(test_args.b, [100])


def test_bad_damping_force(test_args: Args) -> None:
    f_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        damping_def.calculate_velocity(test_args.b, f_bad_vector)

    f_scalar = Quantity(1.0 * units.newton)
    with raises(ValueError):
        damping_def.calculate_velocity(test_args.b, f_scalar)

    with raises(TypeError):
        damping_def.calculate_velocity(test_args.b, 100)
    with raises(TypeError):
        damping_def.calculate_velocity(test_args.b, [100])
