from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, Quantity, errors
from symplyphysics.definitions.vector import momentum_is_mass_times_velocity_vector as momentum_def

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

# Description
## A particle of mass m = 0.2 kg is moving in space, its velocity vector being (-1, 0, 2) m/s.
## The vector of its linear momentum is (-0.2, 0, 0.4) kg*m/s.

Args = namedtuple("Args", "m v p")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(0.2 * units.kilogram)
    v = QuantityCoordinateVector([
        Quantity(-1.0 * units.meter / units.second),
        Quantity(0.0 * units.meter / units.second),
        Quantity(2.0 * units.meter / units.second),
    ], CARTESIAN)
    p = QuantityCoordinateVector([
        Quantity(-0.2 * units.kilogram * units.meter / units.second),
        Quantity(0.0 * units.kilogram * units.meter / units.second),
        Quantity(0.4 * units.kilogram * units.meter / units.second),
    ], CARTESIAN)
    return Args(m=m, v=v, p=p)


def test_momentum_definition(test_args: Args) -> None:
    result = momentum_def.calculate_momentum(test_args.m, test_args.v)
    assert_equal_vectors(result, test_args.p)


def test_momentum_bad_mass(test_args: Args) -> None:
    m_bad_scalar = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(m_bad_scalar, test_args.v)
    with raises(TypeError):
        momentum_def.calculate_momentum(100, test_args.v)


def test_momentum_bad_velocity(test_args: Args) -> None:
    v_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(test_args.m, v_bad_vector)

    v_scalar = Quantity(1.0 * units.meter / units.second)
    with raises(ValueError):
        momentum_def.calculate_momentum(test_args.m, v_scalar)

    with raises(TypeError):
        momentum_def.calculate_momentum(test_args.m, 100)
    with raises(TypeError):
        momentum_def.calculate_momentum(test_args.m, [100, 100])


def test_velocity_law(test_args: Args) -> None:
    result = momentum_def.calculate_velocity(test_args.m, test_args.p)
    assert_equal_vectors(result, test_args.v)


def test_velocity_bad_mass(test_args: Args) -> None:
    m_bad_scalar = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        momentum_def.calculate_velocity(m_bad_scalar, test_args.p)
    with raises(TypeError):
        momentum_def.calculate_velocity(100, test_args.p)


def test_velocity_bad_momentum(test_args: Args) -> None:
    p_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        momentum_def.calculate_velocity(test_args.m, p_bad_vector)

    p_scalar = Quantity(1.0 * units.kilogram * units.meter / units.second)
    with raises(ValueError):
        momentum_def.calculate_velocity(test_args.m, p_scalar)

    with raises(TypeError):
        momentum_def.calculate_velocity(test_args.m, 100)
    with raises(TypeError):
        momentum_def.calculate_velocity(test_args.m, [100])
