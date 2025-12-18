from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.special_relativity.relativistic_kinematics.relativistic_effects import relativistic_length_via_proper_length_and_speed as law

Args = namedtuple("Args", ["l", "v"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    l = Quantity(100 * units.meter)
    v = Quantity(5_000_000 * (units.meter / units.second))
    return Args(l=l, v=v)


def test_basic_length(test_args: Args) -> None:
    result = law.calculate_relativistic_length(test_args.l, test_args.v)
    assert_equal(result, 99.98 * units.meter)


def test_bad_length(test_args: Args) -> None:
    ml = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_relativistic_length(ml, test_args.v)
    with raises(TypeError):
        law.calculate_relativistic_length(100, test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    mv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_relativistic_length(test_args.l, mv)
    with raises(TypeError):
        law.calculate_relativistic_length(test_args.l, 100)
