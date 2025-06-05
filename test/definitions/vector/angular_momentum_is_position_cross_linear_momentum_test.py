from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, Quantity, errors
from symplyphysics.definitions.vector import (
    angular_momentum_is_position_cross_linear_momentum as angular_momentum_def,)

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector, CARTESIAN
from symplyphysics.core.experimental.approx import assert_equal_vectors

# Description
## A particle is located at a point with position vector of (1, 2, -1) m and possesses
## a linear momentum (0, -1, 2) kg*m/s. Its angular momentum amounts to (3, -2, -1) kg*m**2/s.

Args = namedtuple("Args", "r p")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    r = QuantityCoordinateVector([
        Quantity(1.0 * units.meter),
        Quantity(2.0 * units.meter),
        Quantity(-1.0 * units.meter),
    ], CARTESIAN)
    p = QuantityCoordinateVector([
        Quantity(0.0 * units.kilogram * units.meter / units.second),
        Quantity(-1.0 * units.kilogram * units.meter / units.second),
        Quantity(2.0 * units.kilogram * units.meter / units.second),
    ], CARTESIAN)
    return Args(r=r, p=p)


def test_definition(test_args: Args) -> None:
    result = angular_momentum_def.calculate_angular_momentum(test_args.r, test_args.p)
    unit = units.kilogram * units.meter**2 / units.second
    expected = QuantityCoordinateVector([
        3.0 * unit,
        -2.0 * unit,
        -1.0 * unit,
    ], CARTESIAN)
    assert_equal_vectors(result, expected)


def test_bad_position_vector(test_args: Args) -> None:
    r_bad_vector = QuantityCoordinateVector([
        Quantity(0.0 * units.meter / units.second),
        Quantity(-1.0 * units.meter / units.second),
        Quantity(2.0 * units.meter / units.second),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        angular_momentum_def.calculate_angular_momentum(r_bad_vector, test_args.p)

    r_scalar = Quantity(1.0 * units.meter)
    with raises(ValueError):
        angular_momentum_def.calculate_angular_momentum(r_scalar, test_args.p)

    with raises(TypeError):
        angular_momentum_def.calculate_angular_momentum(100, test_args.p)
    with raises(TypeError):
        angular_momentum_def.calculate_angular_momentum([100], test_args.p)


def test_bad_linear_momentum(test_args: Args) -> None:
    p_bad_vector = QuantityCoordinateVector([
        Quantity(0.0 * units.meter / units.second),
        Quantity(-1.0 * units.meter / units.second),
        Quantity(2.0 * units.meter / units.second),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        angular_momentum_def.calculate_angular_momentum(test_args.r, p_bad_vector)

    p_scalar = Quantity(1.0 * units.kilogram * units.meter / units.second)
    with raises(ValueError):
        angular_momentum_def.calculate_angular_momentum(test_args.r, p_scalar)

    with raises(TypeError):
        angular_momentum_def.calculate_angular_momentum(test_args.r, 100)
    with raises(TypeError):
        angular_momentum_def.calculate_angular_momentum(test_args.r, [100])
