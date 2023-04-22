from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import power_is_energy_derivative as power_def


@fixture
def test_args():
    Q0 = Quantity(0 * units.joule)
    Q1 = Quantity(20 * units.joule)
    t = Quantity(5 * units.second)
    Args = namedtuple("Args", ["Q0", "Q1", "t"])
    return Args(Q0=Q0, Q1=Q1, t=t)


def test_basic_power(test_args):
    result = power_def.calculate_power(test_args.Q0, test_args.Q1, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
                                                     units.power)
    result_power = convert_to(result, power_def.definition_units_SI).subs(
        units.watt, 1).evalf(4)
    assert result_power == approx(4, 0.01)


def test_power_with_bad_energy(test_args):
    Qb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        power_def.calculate_power(Qb, test_args.Q1, test_args.t)
    with raises(errors.UnitsError):
        power_def.calculate_power(test_args.Q0, Qb, test_args.t)
    with raises(TypeError):
        power_def.calculate_power(100, test_args.Q1, test_args.t)
    with raises(TypeError):
        power_def.calculate_power(test_args.Q0, 100, test_args.t)


def test_power_with_bad_time(test_args):
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        power_def.calculate_power(test_args.Q0, test_args.Q1, tb)
    with raises(TypeError):
        power_def.calculate_power(test_args.Q0, test_args.Q1, 100)
