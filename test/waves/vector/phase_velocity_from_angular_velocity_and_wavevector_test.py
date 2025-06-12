from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.waves.vector import (
    phase_velocity_from_angular_velocity_and_wavevector as phase_velocity_law,)

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

# Description
## A wave is propagating with the wavevector k = (2, 3, -6) rad/m and angular frequency w = 5 Hz.
## Then its phase speed is (0.204, 0.306, -0.612) m/s.

Args = namedtuple("Args", "v w k")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v = QuantityCoordinateVector([
        Quantity(0.204 * units.meter / units.second),
        Quantity(0.306 * units.meter / units.second),
        Quantity(-0.612 * units.meter / units.second),
    ], CARTESIAN)
    w = Quantity(5.0 * units.hertz)
    k = QuantityCoordinateVector([
        Quantity(2.0 * units.radian / units.meter),
        Quantity(3.0 * units.radian / units.meter),
        Quantity(-6.0 * units.radian / units.meter),
    ], CARTESIAN)
    return Args(v=v, w=w, k=k)


def test_law_velocity(test_args: Args) -> None:
    result = phase_velocity_law.calculate_phase_velocity(test_args.w, test_args.k)
    assert_equal_vectors(result, test_args.v)


def test_law_wavevector(test_args: Args) -> None:
    result = phase_velocity_law.calculate_wavevector(test_args.w, test_args.v)
    assert_equal_vectors(result, test_args.k)


def test_bad_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        phase_velocity_law.calculate_phase_velocity(wb, test_args.k)
    with raises(TypeError):
        phase_velocity_law.calculate_phase_velocity(100, test_args.k)
    with raises(errors.UnitsError):
        phase_velocity_law.calculate_wavevector(wb, test_args.k)
    with raises(TypeError):
        phase_velocity_law.calculate_wavevector(100, test_args.k)


def test_bad_wavenumber(test_args: Args) -> None:
    k_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        phase_velocity_law.calculate_phase_velocity(test_args.w, k_bad_vector)

    k_scalar = Quantity(1.0 * units.radian / units.meter)
    with raises(ValueError):
        phase_velocity_law.calculate_phase_velocity(test_args.w, k_scalar)

    with raises(TypeError):
        phase_velocity_law.calculate_phase_velocity(test_args.w, 100)
    with raises(TypeError):
        phase_velocity_law.calculate_phase_velocity(test_args.w, [100, 100])


def test_bad_velocity(test_args: Args) -> None:
    v_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        phase_velocity_law.calculate_wavevector(test_args.w, v_bad_vector)

    v_scalar = Quantity(1.0 * units.meter / units.second)
    with raises(ValueError):
        phase_velocity_law.calculate_wavevector(test_args.w, v_scalar)

    with raises(TypeError):
        phase_velocity_law.calculate_wavevector(test_args.w, 100)
    with raises(TypeError):
        phase_velocity_law.calculate_wavevector(test_args.w, [100, 100])
