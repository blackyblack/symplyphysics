from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors, S
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.electricity.circuits import sum_of_all_voltages_in_loop_is_zero as kirchhoff_law_2

@fixture
def test_args():
    U1 = Quantity(units.voltage, 3 * units.volt)
    U2 = Quantity(units.voltage, -5 * units.volt)
    Args = namedtuple("Args", ["U1", "U2"])
    return Args(U1=U1, U2=U2)

def test_basic_voltage(test_args):
    result = kirchhoff_law_2.calculate_voltage(test_args.U1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.voltage)
    result_voltage = convert_to(result, units.volt).subs(units.volt, 1).evalf(2)
    assert result_voltage == approx(-3, 0.01)

def test_bad_voltage():
    Ub = Quantity(units.length)
    with raises(errors.UnitsError):
        kirchhoff_law_2.calculate_voltage(Ub)
    with raises(TypeError):
        kirchhoff_law_2.calculate_voltage(100)

def test_array_voltage(test_args):
    result = kirchhoff_law_2.calculate_voltage_from_array([test_args.U1, test_args.U2])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.voltage)
    result_voltage = convert_to(result, units.volt).subs(units.volt, 1).evalf(2)
    assert result_voltage == approx(2, 0.01)

def test_array_empty():
    result = kirchhoff_law_2.calculate_voltage_from_array([])
    assert SI.get_dimension_system().is_dimensionless(result.dimension)
    assert int(convert_to(result, S.One).n()) == 0

def test_array_bad_voltage(test_args):
    Ub = Quantity(units.length)
    with raises(errors.UnitsError):
        kirchhoff_law_2.calculate_voltage_from_array([test_args.U1, Ub])
    with raises(TypeError):
        kirchhoff_law_2.calculate_voltage_from_array([test_args.U1, 100])
    with raises(errors.UnitsError):
        kirchhoff_law_2.calculate_voltage_from_array([Ub, test_args.U2])
    with raises(TypeError):
        kirchhoff_law_2.calculate_voltage_from_array([100, test_args.U2])
    with raises(errors.UnitsError):
        kirchhoff_law_2.calculate_voltage_from_array([Ub, Ub])
    with raises(TypeError):
        kirchhoff_law_2.calculate_voltage_from_array([100, 100])
