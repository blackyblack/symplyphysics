from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, errors, Quantity
from symplyphysics.laws.electricity.vector import (
    force_acting_on_dipole_in_non_uniform_electric_field as law)

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector, CARTESIAN
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "p dx de")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(1e-10 * units.coulomb * units.meter)

    dx = Quantity(1e-3 * units.meter)

    de = QuantityCoordinateVector([
        1 * units.volt / units.meter,
        -3 * units.volt / units.meter,
        2 * units.volt / units.meter,
    ], CARTESIAN)

    return Args(p=p, dx=dx, de=de)


def test_law(test_args: Args) -> None:
    result = law.calculate_force(test_args.p, test_args.dx, test_args.de)

    expected = QuantityCoordinateVector([
        1e-7 * units.newton,
        -3e-7 * units.newton,
        2e-7 * units.newton,
    ], CARTESIAN)

    assert_equal_vectors(result, expected)


def test_bad_electric_dipole_moment(test_args: Args) -> None:
    pb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_force(pb, test_args.dx, test_args.de)
    with raises(TypeError):
        law.calculate_force(100, test_args.dx, test_args.de)


def test_bad_position(test_args: Args) -> None:
    xb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_force(test_args.p, xb, test_args.de)
    with raises(TypeError):
        law.calculate_force(test_args.p, 100, test_args.de)


def test_bad_electric_field(test_args: Args) -> None:
    eb_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_force(test_args.p, test_args.dx, eb_vector)

    eb_scalar = units.volt / units.meter
    with raises(ValueError):
        law.calculate_force(test_args.p, test_args.dx, eb_scalar)

    with raises(TypeError):
        law.calculate_force(test_args.p, test_args.dx, 100)
    with raises(TypeError):
        law.calculate_force(test_args.p, test_args.dx, [100])
