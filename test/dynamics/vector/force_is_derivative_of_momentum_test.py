from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, errors, Quantity
from symplyphysics.laws.dynamics.vector import force_is_derivative_of_momentum as force_momentum_law

from symplyphysics.core.experimental.approx import assert_equal_vectors
from symplyphysics.core.experimental.coordinate_systems import (CartesianCoordinateSystem,
    QuantityCoordinateVector)

# Description
## A force is acting upon a body. At one moment of time, its momentum amounts to (-1, 2, 3) kg*m/s.
## After a second, it was (3, 4, 5) kg*m/s. The force acting on the body amounts to (4, 2, 2) N.

Args = namedtuple("Args", "c p0 p1 dt")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    c = CartesianCoordinateSystem()
    p0 = QuantityCoordinateVector([
        Quantity(-1.0 * units.kilogram * units.meter / units.second),
        Quantity(2.0 * units.kilogram * units.meter / units.second),
        Quantity(3.0 * units.kilogram * units.meter / units.second),
    ], c)
    p1 = QuantityCoordinateVector([
        Quantity(3.0 * units.kilogram * units.meter / units.second),
        Quantity(4.0 * units.kilogram * units.meter / units.second),
        Quantity(5.0 * units.kilogram * units.meter / units.second),
    ], c)
    dt = Quantity(1.0 * units.second)
    return Args(c=c, p0=p0, p1=p1, dt=dt)


def test_basic_law(test_args: Args) -> None:
    result_force = force_momentum_law.calculate_force(test_args.p0, test_args.p1, test_args.dt)

    expected = QuantityCoordinateVector([
        4 * units.newton,
        2 * units.newton,
        2 * units.newton,
    ], test_args.c)

    assert_equal_vectors(result_force, expected)


def test_bad_momenta(test_args: Args) -> None:
    p_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
    ], test_args.c)
    with raises(errors.UnitsError):
        force_momentum_law.calculate_force(p_bad_vector, test_args.p1, test_args.dt)
    with raises(errors.UnitsError):
        force_momentum_law.calculate_force(test_args.p0, p_bad_vector, test_args.dt)

    p_scalar = Quantity(1.0 * units.kilogram * units.meter / units.second)
    with raises(ValueError):
        force_momentum_law.calculate_force(p_scalar, test_args.p1, test_args.dt)
    with raises(ValueError):
        force_momentum_law.calculate_force(test_args.p0, p_scalar, test_args.dt)

    with raises(AttributeError):
        force_momentum_law.calculate_force(
            100,  # type: ignore[arg-type]
            test_args.p1,
            test_args.dt)
    with raises(AttributeError):
        force_momentum_law.calculate_force(
            [100],  # type: ignore[arg-type]
            test_args.p1,
            test_args.dt)
    with raises(AttributeError):
        force_momentum_law.calculate_force(
            test_args.p0,
            100,  # type: ignore[arg-type]
            test_args.dt)
    with raises(AttributeError):
        force_momentum_law.calculate_force(
            test_args.p0,
            [100],  # type: ignore[arg-type]
            test_args.dt)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        force_momentum_law.calculate_force(test_args.p0, test_args.p1, tb)
    with raises(TypeError):
        force_momentum_law.calculate_force(
            test_args.p0,
            test_args.p1,
            100,  # type: ignore[arg-type]
        )
