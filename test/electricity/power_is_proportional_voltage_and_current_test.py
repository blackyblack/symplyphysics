from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import power_is_proportional_voltage_and_current as power_law


# What is the power can be released by the battery at 9 V and a current of 0.5 A?
@fixture
def test_args():
    I = Quantity(0.5 * units.ampere)
    U = Quantity(9 * units.volt)
    Args = namedtuple("Args", ["U", "I"])
    return Args(I=I, U=U)


def test_basic_power(test_args):
    result = power_law.calculate_power(test_args.I, test_args.U)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.power)
    result_current = convert_to(result, units.watt).subs(units.watt, 1).evalf(2)
    assert result_current == approx(4.5, 0.01)


def test_bad_current(test_args):
    Ib = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        power_law.calculate_power(Ib, test_args.U)
    with raises(TypeError):
        power_law.calculate_power(100, test_args.U)


def test_bad_voltage(test_args):
    Ub = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        power_law.calculate_power(test_args.I, Ub)
    with raises(TypeError):
        power_law.calculate_power(test_args.I, 100)
