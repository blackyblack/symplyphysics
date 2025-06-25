from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, errors, Quantity
from symplyphysics.laws.electricity.vector import (
    current_density_is_charge_density_times_drift_velocity as law)

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "rho u")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    rho = Quantity(1e-4 * units.coulomb / units.meter**3)

    u = QuantityCoordinateVector([
        1e-3 * units.meter / units.second,
        0,
        -3e-3 * units.meter / units.second,
    ], CARTESIAN)

    return Args(rho=rho, u=u)


def test_law(test_args: Args) -> None:
    result = law.calculate_current_density(test_args.rho, test_args.u)

    expected = QuantityCoordinateVector([
        1e-7 * units.ampere / units.meter**2,
        0,
        -3e-7 * units.ampere / units.meter**2,
    ], CARTESIAN)

    assert_equal_vectors(result, expected)


def test_bad_charge_density(test_args: Args) -> None:
    rhob = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_current_density(rhob, test_args.u)
    with raises(TypeError):
        law.calculate_current_density(100, test_args.u)


def test_bad_drift_velocity(test_args: Args) -> None:
    ub_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_current_density(test_args.rho, ub_vector)

    ub_scalar = units.meter / units.second
    with raises(ValueError):
        law.calculate_current_density(test_args.rho, ub_scalar)

    with raises(TypeError):
        law.calculate_current_density(test_args.rho, 100)
    with raises(TypeError):
        law.calculate_current_density(test_args.rho, [100])
