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
from symplyphysics.laws.electricity import power_is_proportional_voltage_and_current as power_law


# What is the power can be released by the battery at 9 V and a current of 0.5 A?
@fixture(name="test_args")
def test_args_fixture():
    I = Quantity(0.5 * units.ampere)
    U = Quantity(9 * units.volt)
    Args = namedtuple("Args", ["U", "I"])
    return Args(I=I, U=U)


def test_basic_power(test_args):
    result = power_law.calculate_power(test_args.I, test_args.U)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.power)
    result_power = convert_to(result, units.watt).evalf(2)
    assert_approx(result_power, 4.5)


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
