from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import S
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity.circuits import sum_of_all_voltages_in_loop_is_zero as kirchhoff_law_2


@fixture(name="test_args")
def test_args_fixture():
    U1 = Quantity(3 * units.volt)
    U2 = Quantity(-5 * units.volt)
    Args = namedtuple("Args", ["U1", "U2"])
    return Args(U1=U1, U2=U2)


def test_basic_voltage(test_args):
    result = kirchhoff_law_2.calculate_voltage([test_args.U1])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.voltage)
    result_voltage = convert_to(result, units.volt).evalf(2)
    assert result_voltage == approx(-3, 0.01)


def test_three_voltage_array(test_args):
    result = kirchhoff_law_2.calculate_voltage([test_args.U1, test_args.U2])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.voltage)
    result_voltage = convert_to(result, units.volt).evalf(2)
    assert result_voltage == approx(2, 0.01)


def test_array_empty():
    result = kirchhoff_law_2.calculate_voltage([])
    assert SI.get_dimension_system().is_dimensionless(result.dimension)
    assert int(convert_to(result, S.One).n()) == 0


def test_array_bad_voltage(test_args):
    Ub = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        kirchhoff_law_2.calculate_voltage([test_args.U1, Ub])
    with raises(TypeError):
        kirchhoff_law_2.calculate_voltage([test_args.U1, 100])
    with raises(errors.UnitsError):
        kirchhoff_law_2.calculate_voltage([Ub, test_args.U2])
    with raises(TypeError):
        kirchhoff_law_2.calculate_voltage([100, test_args.U2])
    with raises(errors.UnitsError):
        kirchhoff_law_2.calculate_voltage([Ub, Ub])
    with raises(TypeError):
        kirchhoff_law_2.calculate_voltage([100, 100])
