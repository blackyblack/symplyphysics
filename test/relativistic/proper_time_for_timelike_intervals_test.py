from collections import namedtuple
from pytest import fixture, raises
from sympy import I
from symplyphysics import (
    assert_equal,
    errors,
    prefixes,
    units,
    Quantity,
)
from symplyphysics.laws.relativistic import proper_time_for_timelike_intervals as time_law

Args = namedtuple("Args", "ds")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    ds = Quantity(1 * units.millimeter)
    return Args(ds=ds)


def test_law(test_args: Args) -> None:
    result = time_law.calculate_proper_time(test_args.ds)
    assert_equal(result, 3.34 * prefixes.pico * units.second, tolerance=2e-3)


def test_bad_interval() -> None:
    sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        time_law.calculate_proper_time(sb)
    with raises(TypeError):
        time_law.calculate_proper_time(100)

    sb = Quantity((1 + 2 * I) * units.meter)
    with raises(ValueError):
        time_law.calculate_proper_time(sb)
