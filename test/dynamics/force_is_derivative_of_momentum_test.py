from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, units, Quantity, errors
from symplyphysics.laws.dynamics import force_is_derivative_of_momentum as law

Args = namedtuple("Args", "dp dt")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dp = Quantity(6 * units.kilogram * units.meter / units.second)
    dt = Quantity(2 * units.second)
    return Args(dp=dp, dt=dt)


def test_law(test_args: Args) -> None:
    result = law.calculate_force(test_args.dp, test_args.dt)
    assert_equal(result, 3 * units.newton)


def test_bad_momentum(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_force(pb, test_args.dt)
    with raises(TypeError):
        law.calculate_force(100, test_args.dt)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_force(test_args.dp, tb)
    with raises(TypeError):
        law.calculate_force(test_args.dp, 100)
