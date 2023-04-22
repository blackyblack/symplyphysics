from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import period_from_circular_frequency as period_def


@fixture
def test_args():
    w = Quantity(6.28 * units.radian / units.second)
    Args = namedtuple("Args", ["w"])
    return Args(w=w)


def test_basic_period(test_args):
    result = period_def.calculate_period(test_args.w)
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
                                                     units.time)
    result_period = convert_to(result, period_def.definition_units_SI).subs(
        units.second, 1).evalf(2)
    assert result_period == approx(1.0, 0.01)


def test_bad_frequency():
    wb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        period_def.calculate_period(wb)
    with raises(TypeError):
        period_def.calculate_period(100)
