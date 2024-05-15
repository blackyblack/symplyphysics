from collections import namedtuple
from pytest import fixture, raises
from sympy import I
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.relativistic import spacetime_interval as law

Args = namedtuple("Args", "dt dx")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dt = Quantity(1 * units.nanosecond)
    dx = Quantity(1 * units.kilometer)
    return Args(dt=dt, dx=dx)


def test_law(test_args: Args) -> None:
    result = law.calculate_spacetime_interval(test_args.dt, test_args.dx)
    assert_equal(result, 1000 * I * units.meter)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_spacetime_interval(tb, test_args.dx)
    with raises(TypeError):
        law.calculate_spacetime_interval(100, test_args.dx)


def test_bad_position(test_args: Args) -> None:
    xb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_spacetime_interval(test_args.dt, xb)
    with raises(TypeError):
        law.calculate_spacetime_interval(test_args.dt, 100)
