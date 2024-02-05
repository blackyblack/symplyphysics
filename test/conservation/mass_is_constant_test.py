from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.conservation import mass_is_constant as conservation_law


@fixture(name="test_args")
def test_args_fixture():
    ms = Quantity(5 * units.kilograms)
    Args = namedtuple("Args", ["ms"])
    return Args(ms=ms)


def test_basic_conservation(test_args):
    result_expr = conservation_law.calculate_mass_after(test_args.ms)
    assert SI.get_dimension_system().equivalent_dims(result_expr.dimension, units.mass)
    result = convert_to(result_expr, units.kilograms).evalf(2)
    assert_approx(result, 5)


def test_bad_mass():
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        conservation_law.calculate_mass_after(mb)
    with raises(TypeError):
        conservation_law.calculate_mass_after(100)
