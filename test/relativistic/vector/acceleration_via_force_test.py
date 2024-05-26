from collections import namedtuple
from pytest import fixture, raises
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Quantity,
    QuantityVector,
    errors,
)
from symplyphysics.core.approx import assert_equal_vectors
from symplyphysics.laws.relativistic.vector import acceleration_via_force as law

Args = namedtuple("Args", "m f v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(units.electron_rest_mass)
    f = QuantityVector([
        Quantity(0),
        Quantity(-4.13e-27 * units.newton),
        Quantity(2.75e-27 * units.newton),
    ])
    v = QuantityVector([
        Quantity(speed_of_light / 2),
        Quantity(speed_of_light / 2),
        Quantity(-1 * speed_of_light / 4),
    ])
    return Args(m=m, f=f, v=v)


def test_law(test_args: Args) -> None:
    result = law.calculate_acceleration(test_args.m, test_args.f, test_args.v)
    correct = QuantityVector([
        Quantity(1.0 * units.kilometer / units.second**2),
        Quantity(-2.0 * units.kilometer / units.second**2),
        Quantity(1.5 * units.kilometer / units.second**2),
    ])
    assert_equal_vectors(result, correct, tolerance=2e-3)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_acceleration(mb, test_args.f, test_args.v)
    with raises(TypeError):
        law.calculate_acceleration(100, test_args.f, test_args.v)


def test_bad_force(test_args: Args) -> None:
    fb_vector = QuantityVector([Quantity(units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.m, fb_vector, test_args.v)
    
    fb_scalar = Quantity(units.newton)
    with raises(AttributeError):
        law.calculate_acceleration(test_args.m, fb_scalar, test_args.v)

    with raises(TypeError):
        law.calculate_acceleration(test_args.m, 100, test_args.v)
    with raises(TypeError):
        law.calculate_acceleration(test_args.m, [100], test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityVector([Quantity(units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.m, test_args.f, vb_vector)
    
    vb_scalar = Quantity(units.meter / units.second)
    with raises(AttributeError):
        law.calculate_acceleration(test_args.m, test_args.f, vb_scalar)

    with raises(TypeError):
        law.calculate_acceleration(test_args.m, test_args.f, 100)
    with raises(TypeError):
        law.calculate_acceleration(test_args.m, test_args.f, [100])
