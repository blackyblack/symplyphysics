from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import momentum_derivative_of_kinetic_energy_is_speed as law

Args = namedtuple("Args", "de dp")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    de = Quantity(1e-3 * units.joule)
    dp = Quantity(2e-4 * units.kilogram * units.meter / units.second)
    return Args(de=de, dp=dp)


def test_law(test_args: Args) -> None:
    result = law.calculate_speed(test_args.de, test_args.dp)
    assert_equal(result, 5 * units.meter / units.second)


def test_bad_energy(test_args: Args) -> None:
    eb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_speed(eb, test_args.dp)
    with raises(TypeError):
        law.calculate_speed(100, test_args.dp)


def test_bad_momentum(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_speed(test_args.de, pb)
    with raises(TypeError):
        law.calculate_speed(test_args.de, 100)
