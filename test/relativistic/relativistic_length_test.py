from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    convert_to,
)
from symplyphysics.laws.relativistic import relativistic_length


@fixture(name="test_args")
def test_args_fixture():
    l = Quantity(100 * units.meters)
    v = Quantity(5_000_000 * (units.meter / units.second))
    Args = namedtuple("Args", ["l", "v"])
    return Args(l=l, v=v)


def test_basic_length(test_args):
    result = relativistic_length.calculate_relativistic_length(
        test_args.l, test_args.v)
    result_length = convert_to(result, units.length).evalf(4)
    assert result_length == approx(99.98, 0.001)


def test_bad_length(test_args):
    ml = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relativistic_length.calculate_relativistic_length(ml, test_args.v)
    with raises(TypeError):
        relativistic_length.calculate_relativistic_length(100, test_args.v)


def test_bad_velocity(test_args):
    mv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relativistic_length.calculate_relativistic_length(test_args.l, mv)
    with raises(TypeError):
        relativistic_length.calculate_relativistic_length(test_args.l, 100)
