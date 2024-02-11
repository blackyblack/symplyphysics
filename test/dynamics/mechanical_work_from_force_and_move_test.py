from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import mechanical_work_from_force_and_move as work_law

# Description
## Force of 100N is applied to heavy object lying oh horizontal table.
## Object slides on the table with no friction and is been moved to 3m far.
## Work should be 100 * 3 = 300 Joules.

Args = namedtuple("Args", ["F", "S"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    F = Quantity(100 * units.newton)
    S = Quantity(3 * units.meter)
    return Args(F=F, S=S)


def test_basic_work(test_args: Args) -> None:
    result = work_law.calculate_work(test_args.F, test_args.S)
    assert_equal(result, 300 * units.joule)


def test_bad_force(test_args: Args) -> None:
    Fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        work_law.calculate_work(Fb, test_args.S)
    with raises(TypeError):
        work_law.calculate_work(100, test_args.S)


def test_bad_move(test_args: Args) -> None:
    Sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.F, Sb)
    with raises(TypeError):
        work_law.calculate_work(test_args.F, 100)
