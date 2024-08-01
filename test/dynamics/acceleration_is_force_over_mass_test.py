from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass as newton_second_law

Args = namedtuple("Args", ["m", "a"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(1 * units.kilogram)
    a = Quantity(3 * units.meter / units.second**2)
    return Args(m=m, a=a)


def test_basic_force(test_args: Args) -> None:
    result = newton_second_law.calculate_force(test_args.m, test_args.a)
    assert_equal(result, 3 * units.newton)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_force(mb, test_args.a)
    with raises(TypeError):
        newton_second_law.calculate_force(100, test_args.a)


def test_bad_acceleration(test_args: Args) -> None:
    ab = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_second_law.calculate_force(test_args.m, ab)
    with raises(TypeError):
        newton_second_law.calculate_force(test_args.m, 100)
