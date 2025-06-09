from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, errors, Quantity
from symplyphysics.laws.dynamics.vector import (torque_is_angular_momentum_derivative as torque_law)

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

# Description
## During the interval of 1 s, the particle's angular momentum changed from (1, 2, -1) kg*m**2/s
## to (0, 2, 1) kg*m**2/s. The torque acting on the particle was (-1, 0, 2) N*m.

Args = namedtuple("Args", "t L0 L1")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(1.0 * units.second)
    L0 = QuantityCoordinateVector([
        Quantity(1.0 * units.kilogram * units.meter**2 / units.second),
        Quantity(2.0 * units.kilogram * units.meter**2 / units.second),
        Quantity(-1.0 * units.kilogram * units.meter**2 / units.second),
    ], CARTESIAN)
    L1 = QuantityCoordinateVector([
        Quantity(0.0 * units.kilogram * units.meter**2 / units.second),
        Quantity(2.0 * units.kilogram * units.meter**2 / units.second),
        Quantity(1.0 * units.kilogram * units.meter**2 / units.second),
    ], CARTESIAN)
    return Args(t=t, L0=L0, L1=L1)


def test_law(test_args: Args) -> None:
    result = torque_law.calculate_torque(test_args.L0, test_args.L1, test_args.t)
    expected = QuantityCoordinateVector([
        -1 * units.newton * units.meter,
        0,
        2 * units.newton * units.meter,
    ], CARTESIAN)
    assert_equal_vectors(result, expected)


def test_bad_angular_momenta(test_args: Args) -> None:
    Lb = QuantityCoordinateVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        torque_law.calculate_torque(Lb, test_args.L1, test_args.t)
    with raises(errors.UnitsError):
        torque_law.calculate_torque(test_args.L0, Lb, test_args.t)

    Ls = Quantity(1.0 * units.kilogram * units.meter**2 / units.second)
    with raises(TypeError):
        torque_law.calculate_torque(Ls, test_args.L1, test_args.t)
    with raises(TypeError):
        torque_law.calculate_torque(test_args.L0, Ls, test_args.t)

    with raises(TypeError):
        torque_law.calculate_torque(100, test_args.L1, test_args.t)
    with raises(TypeError):
        torque_law.calculate_torque([100], test_args.L1, test_args.t)
    with raises(TypeError):
        torque_law.calculate_torque(test_args.L0, 100, test_args.t)
    with raises(TypeError):
        torque_law.calculate_torque(test_args.L0, [100], test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        torque_law.calculate_torque(test_args.L0, test_args.L1, tb)
    with raises(TypeError):
        torque_law.calculate_torque(test_args.L0, test_args.L1, 100)
