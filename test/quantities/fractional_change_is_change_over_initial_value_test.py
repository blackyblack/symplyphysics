from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.quantities import (
    fractional_change_is_change_over_initial_value as law,
)

Args = namedtuple("Args", "dm m dv v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dm = Quantity(0.010 * units.kilogram)
    m = Quantity(5.0 * units.kilogram)
    dv = Quantity(0.001 * units.liter)
    v = Quantity(2 * units.liter)
    return Args(dm=dm, m=m, dv=dv, v=v)


def test_mass_law(test_args: Args) -> None:
    result = law.calculate_fractional_change(test_args.dm, test_args.m)
    assert_equal(result, 2.0e-3)


def test_volume_law(test_args: Args) -> None:
    result = law.calculate_fractional_change(test_args.dv, test_args.v)
    assert_equal(result, 5e-4)


def test_bad_arguments(test_args: Args) -> None:
    with raises(errors.UnitsError):
        law.calculate_fractional_change(test_args.dm, test_args.v)
    with raises(errors.UnitsError):
        law.calculate_fractional_change(test_args.dv, test_args.m)
