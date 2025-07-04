from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.dynamics.vector import (restoring_torque_due_to_twist_of_torsion_pendulum as
    pendulum_laws)

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector, CARTESIAN
from symplyphysics.core.experimental.approx import assert_equal_vectors

# Description
## A torsion pendulum is rotated about its axis so that its rotation pseudovector is (0.1, 0.2, -0.1) rad.
## The torsion constant of the pendulum is 10.0 N*m. Vector of the restoring torque of the pendulum amounts
## to (-1.0, -2.0, 1.0) N*m.

Args = namedtuple("Args", "tau kappa theta")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    tau = QuantityCoordinateVector([
        Quantity(-1.0 * units.newton * units.meter),
        Quantity(-2.0 * units.newton * units.meter),
        Quantity(1.0 * units.newton * units.meter),
    ], CARTESIAN)
    kappa = Quantity(10.0 * units.newton * units.meter)
    theta = QuantityCoordinateVector([0.1, 0.2, -0.1], CARTESIAN)
    return Args(tau=tau, kappa=kappa, theta=theta)


def test_torque_law(test_args: Args) -> None:
    result = pendulum_laws.calculate_torque(test_args.kappa, test_args.theta)
    assert_equal_vectors(result, test_args.tau)


def test_rotation_vector_law(test_args: Args) -> None:
    result = pendulum_laws.calculate_rotation_vector(test_args.kappa, test_args.tau)
    assert_equal_vectors(result, test_args.theta)


def test_bad_torsion_constant(test_args: Args) -> None:
    kappab = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        pendulum_laws.calculate_torque(kappab, test_args.theta)
    with raises(TypeError):
        pendulum_laws.calculate_torque(100, test_args.theta)
    with raises(errors.UnitsError):
        pendulum_laws.calculate_rotation_vector(kappab, test_args.tau)
    with raises(TypeError):
        pendulum_laws.calculate_rotation_vector(100, test_args.tau)


def test_bad_rotation_vector(test_args: Args) -> None:
    theta_b = QuantityCoordinateVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        pendulum_laws.calculate_torque(test_args.kappa, theta_b)

    theta_s = 1.0
    with raises(ValueError):
        pendulum_laws.calculate_torque(test_args.kappa, theta_s)


def test_bad_torque(test_args: Args) -> None:
    tau_b = QuantityCoordinateVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        pendulum_laws.calculate_rotation_vector(test_args.kappa, tau_b)

    tau_s = Quantity(1.0 * units.newton * units.meter)
    with raises(ValueError):
        pendulum_laws.calculate_rotation_vector(test_args.kappa, tau_s)

    with raises(TypeError):
        pendulum_laws.calculate_rotation_vector(test_args.kappa, 100)
    with raises(TypeError):
        pendulum_laws.calculate_rotation_vector(test_args.kappa, [100])
