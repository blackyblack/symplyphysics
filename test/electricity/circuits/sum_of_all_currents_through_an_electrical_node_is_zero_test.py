from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import S
from symplyphysics import (
    errors, units, convert_to, Quantity, SI,
)
from symplyphysics.laws.electricity.circuits import sum_of_all_currents_through_an_electrical_node_is_zero as kirchhoff_law

@fixture
def test_args():
    I1 = Quantity(3 * units.ampere)
    I2 = Quantity(-5 * units.ampere)
    Args = namedtuple("Args", ["I1", "I2"])
    return Args(I1=I1, I2=I2)

def test_basic_current(test_args):
    result = kirchhoff_law.calculate_current(test_args.I1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.current)
    result_current = convert_to(result, units.ampere).subs(units.ampere, 1).evalf(2)
    assert result_current == approx(-3, 0.01)

def test_bad_current():
    Ib = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        kirchhoff_law.calculate_current(Ib)
    with raises(TypeError):
        kirchhoff_law.calculate_current(100)

def test_array_current(test_args):
    result = kirchhoff_law.calculate_current_from_array([test_args.I1, test_args.I2])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.current)
    result_current = convert_to(result, units.ampere).subs(units.ampere, 1).evalf(2)
    assert result_current == approx(2, 0.01)

def test_array_empty():
    result = kirchhoff_law.calculate_current_from_array([])
    assert SI.get_dimension_system().is_dimensionless(result.dimension)
    assert int(convert_to(result, S.One).n()) == 0

def test_array_bad_current(test_args):
    Ib = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        kirchhoff_law.calculate_current_from_array([test_args.I1, Ib])
    with raises(TypeError):
        kirchhoff_law.calculate_current_from_array([test_args.I1, 100])
    with raises(errors.UnitsError):
        kirchhoff_law.calculate_current_from_array([Ib, test_args.I2])
    with raises(TypeError):
        kirchhoff_law.calculate_current_from_array([100, test_args.I2])
    with raises(errors.UnitsError):
        kirchhoff_law.calculate_current_from_array([Ib, Ib])
    with raises(TypeError):
        kirchhoff_law.calculate_current_from_array([100, 100])
