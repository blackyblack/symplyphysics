from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import moment_of_force_from_moment_of_inertia_and_angular_acceleration as moment_force_law

Args = namedtuple("Args", ["i", "epsilon"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    epsilon = Quantity(1 * units.radians / (units.second**2))
    i = Quantity(3 * units.kilograms * units.meters**2)
    return Args(i=i, epsilon=epsilon)


def test_basic_moment_of_force(test_args: Args) -> None:
    result = moment_force_law.calculate_moment_of_force(test_args.i, test_args.epsilon)
    assert_equal(result, 3 * units.newtons * units.meters)


def test_bad_moment_of_inertia(test_args: Args) -> None:
    ib = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        moment_force_law.calculate_moment_of_force(ib, test_args.epsilon)
    with raises(TypeError):
        moment_force_law.calculate_moment_of_force(100, test_args.i)


def test_bad_angular_acceleration(test_args: Args) -> None:
    epsilon_b = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        moment_force_law.calculate_moment_of_force(test_args.i, epsilon_b)
    with raises(TypeError):
        moment_force_law.calculate_moment_of_force(test_args.i, 100)
