from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, Quantity, errors
from symplyphysics.quantities import speed_of_light
from symplyphysics.laws.relativistic.vector import acceleration_from_force_and_velocity as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "m a v f")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(units.electron_rest_mass)
    a = QuantityCoordinateVector([
        Quantity(1.0 * units.kilometer / units.second**2),
        Quantity(-2.0 * units.kilometer / units.second**2),
        Quantity(1.5 * units.kilometer / units.second**2),
    ], CARTESIAN)
    v = QuantityCoordinateVector([
        Quantity(speed_of_light / 2),
        Quantity(speed_of_light / 2),
        Quantity(-1 * speed_of_light / 4),
    ], CARTESIAN)
    f = QuantityCoordinateVector([
        Quantity(0),
        Quantity(-4.13e-27 * units.newton),
        Quantity(2.75e-27 * units.newton),
    ], CARTESIAN)
    return Args(m=m, a=a, v=v, f=f)


def test_law(test_args: Args) -> None:
    result = law.calculate_acceleration(test_args.m, test_args.f, test_args.v)
    assert_equal_vectors(result, test_args.a, relative_tolerance=2e-3)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_acceleration(mb, test_args.f, test_args.v)
    with raises(TypeError):
        law.calculate_acceleration(100, test_args.f, test_args.v)


def test_bad_force(test_args: Args) -> None:
    fb_vector = QuantityCoordinateVector([Quantity(units.coulomb), 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.m, fb_vector, test_args.v)

    fb_scalar = Quantity(units.planck_force)
    with raises(ValueError):
        law.calculate_acceleration(test_args.m, fb_scalar, test_args.v)

    with raises(TypeError):
        law.calculate_acceleration(test_args.m, 100, test_args.v)
    with raises(TypeError):
        law.calculate_acceleration(test_args.m, [100], test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityCoordinateVector([Quantity(units.coulomb), 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.m, test_args.f, vb_vector)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(ValueError):
        law.calculate_acceleration(test_args.m, test_args.f, vb_scalar)

    with raises(TypeError):
        law.calculate_acceleration(test_args.m, test_args.f, 100)
    with raises(TypeError):
        law.calculate_acceleration(test_args.m, test_args.f, [100])
