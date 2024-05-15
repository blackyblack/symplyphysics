from collections import namedtuple
from pytest import fixture, raises
from sympy import I
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.relativistic import spacetime_interval_is_lorentz_invariant as invariant_law

Args = namedtuple("Args", "s")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    s = Quantity(5 * I * units.meter)
    return Args(s=s)


def test_law(test_args: Args) -> None:
    result = invariant_law.calculate_second_spacetime_interval(test_args.s)
    assert_equal(result, 5 * I * units.meter)


def test_bad_interval() -> None:
    sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        invariant_law.calculate_second_spacetime_interval(sb)
    with raises(TypeError):
        invariant_law.calculate_second_spacetime_interval(100)
