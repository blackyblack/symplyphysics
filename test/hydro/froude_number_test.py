from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, Quantity, units, errors
from symplyphysics.laws.hydro import froude_number

# Example from https://www.symbolab.com/calculator/physics/froude-number

Args = namedtuple("Args", ["l", "v"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    l = Quantity(1 * units.meter)
    v = Quantity(0.5 * units.meter / units.second)
    return Args(l=l, v=v)


def test_basic_froude_number(test_args: Args) -> None:
    result = froude_number.calculate_froude_number(test_args.v, test_args.l)
    assert_equal(result, 0.1596)


def test_bad_velocity(test_args: Args) -> None:
    bd = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        froude_number.calculate_froude_number(bd, test_args.l)
    with raises(TypeError):
        froude_number.calculate_froude_number(10, test_args.l)


def test_bad_length(test_args: Args) -> None:
    bd = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        froude_number.calculate_froude_number(test_args.v, bd)
    with raises(TypeError):
        froude_number.calculate_froude_number(test_args.v, 10)
