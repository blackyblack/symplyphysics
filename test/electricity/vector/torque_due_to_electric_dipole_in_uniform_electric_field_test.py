from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, errors
from symplyphysics.laws.electricity.vector import (
    torque_due_to_electric_dipole_in_uniform_electric_field as law)

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "p e")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = QuantityCoordinateVector([
        1e-10 * units.coulomb * units.meter,
        0,
        -4e-10 * units.coulomb * units.meter,
    ], CARTESIAN)

    e = QuantityCoordinateVector([
        100 * units.volt / units.meter,
        -90 * units.volt / units.meter,
        -1 * units.volt / units.meter,
    ], CARTESIAN)

    return Args(p=p, e=e)


def test_law(test_args: Args) -> None:
    result = law.calculate_torque(test_args.p, test_args.e)

    expected = QuantityCoordinateVector([
        -3.6e-8 * units.newton * units.meter,
        -4.0e-8 * units.newton * units.meter,
        -9.0e-9 * units.newton * units.meter,
    ], CARTESIAN)

    assert_equal_vectors(result, expected, relative_tolerance=3e-3)


def test_bad_electric_dipole_moment(test_args: Args) -> None:
    pb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_torque(pb_vector, test_args.e)

    pb_scalar = units.meter * units.coulomb
    with raises(ValueError):
        law.calculate_torque(pb_scalar, test_args.e)

    with raises(TypeError):
        law.calculate_torque(100, test_args.e)
    with raises(TypeError):
        law.calculate_torque([100], test_args.e)


def test_bad_electric_field(test_args: Args) -> None:
    eb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_torque(test_args.p, eb_vector)

    eb_scalar = units.volt / units.meter
    with raises(ValueError):
        law.calculate_torque(test_args.p, eb_scalar)

    with raises(TypeError):
        law.calculate_torque(test_args.p, 100)
    with raises(TypeError):
        law.calculate_torque(test_args.p, [100])
