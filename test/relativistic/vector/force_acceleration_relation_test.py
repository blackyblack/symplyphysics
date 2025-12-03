from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, Quantity, errors
from symplyphysics.quantities import speed_of_light
from symplyphysics.laws.relativistic.vector import force_acceleration_relation as law

from symplyphysics.core.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.approx import assert_equal_vectors

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
    result = law.calculate_force(test_args.m, test_args.a, test_args.v)
    assert_equal_vectors(result, test_args.f, absolute_tolerance=1e-20)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_force(mb, test_args.a, test_args.v)
    with raises(TypeError):
        law.calculate_force(100, test_args.a, test_args.v)


def test_bad_acceleration(test_args: Args) -> None:
    ab_vector = QuantityCoordinateVector([Quantity(units.coulomb), 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_force(test_args.m, ab_vector, test_args.v)

    ab_scalar = Quantity(units.planck_acceleration)
    with raises(ValueError):
        law.calculate_force(test_args.m, ab_scalar, test_args.v)

    with raises(TypeError):
        law.calculate_force(test_args.m, 100, test_args.v)
    with raises(TypeError):
        law.calculate_force(test_args.m, [100], test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityCoordinateVector([Quantity(units.coulomb), 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_force(test_args.m, test_args.a, vb_vector)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(ValueError):
        law.calculate_force(test_args.m, test_args.a, vb_scalar)

    with raises(TypeError):
        law.calculate_force(test_args.m, test_args.a, 100)
    with raises(TypeError):
        law.calculate_force(test_args.m, test_args.a, [100])
