from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, errors, Quantity
from symplyphysics.laws.dynamics.vector import acceleration_from_force as newton_second_law

from symplyphysics.core.experimental.coordinate_systems import (CartesianCoordinateSystem,
    QuantityCoordinateVector)
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", ["cart", "m", "a", "f"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cart = CartesianCoordinateSystem()
    m = Quantity(1 * units.kilogram)
    a = QuantityCoordinateVector([Quantity(3 * units.meter / units.second**2), 0, 0], cart)
    f = QuantityCoordinateVector([Quantity(3 * units.newton), 0, 0], cart)
    return Args(cart=cart, m=m, a=a, f=f)


def test_basic_force(test_args: Args) -> None:
    result = newton_second_law.calculate_force(test_args.m, test_args.a)
    expected = QuantityCoordinateVector([3 * units.newton, 0, 0], test_args.cart)
    assert_equal_vectors(result, expected)


def test_basic_acceleration(test_args: Args) -> None:
    result = newton_second_law.calculate_acceleration(test_args.m, test_args.f)
    expected = QuantityCoordinateVector([3 * units.meter / units.second**2, 0, 0], test_args.cart)
    assert_equal_vectors(result, expected)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_force(mb, test_args.a)
    with raises(TypeError):
        newton_second_law.calculate_force(100, test_args.a)
    with raises(errors.UnitsError):
        newton_second_law.calculate_acceleration(mb, test_args.f)
    with raises(TypeError):
        newton_second_law.calculate_acceleration(100, test_args.f)


def test_bad_acceleration(test_args: Args) -> None:
    ab = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_force(test_args.m, ab)
    with raises(TypeError):
        newton_second_law.calculate_force(test_args.m, 100)


def test_bad_force(test_args: Args) -> None:
    fb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_acceleration(test_args.m, fb)
    with raises(TypeError):
        newton_second_law.calculate_acceleration(test_args.m, 100)
