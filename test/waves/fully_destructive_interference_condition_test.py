from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import fully_destructive_interference_condition as condition_law

# Description
## The interference of two waves is fully destructive and the phase shift between them is
## -pi, i.e. n = -1.

Args = namedtuple("Args", "n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n = -1
    return Args(n=n)


def test_law(test_args: Args) -> None:
    result = condition_law.calculate_phase_shift(test_args.n)
    assert_equal(result, -1 * pi)


def test_bad_integer_factor() -> None:
    nb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        condition_law.calculate_phase_shift(nb)
