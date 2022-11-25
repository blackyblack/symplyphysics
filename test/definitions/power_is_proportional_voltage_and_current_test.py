from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.definitions import power_is_proportional_voltage_and_current as power_def
# What is the power can release the battery at 9 V and a current of 0.5 A?
@fixture
def test_args():
    I = units.Quantity('I')
    SI.set_quantity_dimension(I, units.current)
    SI.set_quantity_scale_factor(I, 0.5 * units.ampere)
    U = units.Quantity('U')
    SI.set_quantity_dimension(U, units.voltage)
    SI.set_quantity_scale_factor(U, 9 * units.volt)
    Args = namedtuple('Args', ['U', 'I'])
    return Args(I=I, U=U)

def test_basic_power(test_args):
    result = power_def.calculate_power(test_args.I, test_args.U)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.power)
    result_current = convert_to(result, units.watt).subs(units.watt, 1).evalf(2)
    assert result_current == approx(4.5, 0.01)

def test_bad_current(test_args):
    bI = units.Quantity('bI')
    SI.set_quantity_dimension(bI, units.length)
    SI.set_quantity_scale_factor(bI, 1 * units.meter)

    with raises(errors.UnitsError):
        power_def.calculate_power(bI, test_args.U)

    with raises(TypeError):
        power_def.calculate_power(100, test_args.U)

def test_bad_voltage(test_args):
    bU = units.Quantity('bU')
    SI.set_quantity_dimension(bU, units.length)
    SI.set_quantity_scale_factor(bU, 1 * units.meter)

    with raises(errors.UnitsError):
        power_def.calculate_power(test_args.I, bU)

    with raises(TypeError):
        power_def.calculate_power(test_args.I, 100)