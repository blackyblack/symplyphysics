from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import kinetic_energy_differential as law

Args = namedtuple("Args", "v dp")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v = Quantity(5 * units.meter / units.second)
    dp = Quantity(2e-4 * units.kilogram * units.meter / units.second)
    return Args(v=v, dp=dp)


def test_law(test_args: Args) -> None:
    result = law.calculate_kinetic_energy_differential(test_args.v, test_args.dp)
    assert_equal(result, 1e-3 * units.joule)


def test_bad_speed(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_kinetic_energy_differential(vb, test_args.dp)
    with raises(TypeError):
        law.calculate_kinetic_energy_differential(100, test_args.dp)


def test_bad_momentum(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_kinetic_energy_differential(test_args.v, pb)
    with raises(TypeError):
        law.calculate_kinetic_energy_differential(test_args.v, 100)
