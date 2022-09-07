from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.definitions import current_is_charge_derivative as current_def

@fixture
def test_args():
    Q0 = units.Quantity('Q0')
    SI.set_quantity_dimension(Q0, units.charge)
    SI.set_quantity_scale_factor(Q0, 0 * units.coulomb)
    Q1 = units.Quantity('Q1')
    SI.set_quantity_dimension(Q1, units.charge)
    SI.set_quantity_scale_factor(Q1, 20 * units.coulomb)
    
    I0 = units.Quantity('I0')
    SI.set_quantity_dimension(I0, units.current)
    SI.set_quantity_scale_factor(I0, 0.5 * units.ampere)
    I1 = units.Quantity('I1')
    SI.set_quantity_dimension(I1, units.current)
    SI.set_quantity_scale_factor(I1, 0.6 * units.ampere)

    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.time)
    SI.set_quantity_scale_factor(t, 5 * units.second)

    Args = namedtuple('Args', ['Q0', 'Q1', 'I0', 'I1', 't'])
    return Args(Q0=Q0, Q1=Q1, I0=I0, I1=I1, t=t)


def test_basic_current(test_args):
    result = current_def.calculate_current(
        test_args.Q0, test_args.Q1, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(
        result.dimension, units.current)

    result_current = convert_to(result, current_def.definition_dimension_SI).subs({
        units.ampere: 1}).evalf(2)
    assert result_current == approx(4, 0.01)


def test_current_with_bad_charge(test_args):
    Q0b = units.Quantity('Q0b')
    SI.set_quantity_dimension(Q0b, units.length)
    SI.set_quantity_scale_factor(Q0b, 1 * units.meter)
    with raises(errors.UnitsError):
        current_def.calculate_current(
            Q0b, test_args.Q1, test_args.t)

    # Let Q1 be invalid    
    with raises(errors.UnitsError):
        current_def.calculate_current(
            test_args.Q0, Q0b, test_args.t)

    with raises(TypeError):
        current_def.calculate_current(
            100, test_args.Q1, test_args.t)

    with raises(TypeError):
        current_def.calculate_current(
            test_args.Q0, 100, test_args.t)


def test_current_with_bad_time(test_args):
    tb = units.Quantity('tb')
    SI.set_quantity_dimension(tb, units.length)
    SI.set_quantity_scale_factor(tb, 1 * units.meter)

    with raises(errors.UnitsError):
        current_def.calculate_current(
            test_args.Q0, test_args.Q1, tb)

    with raises(TypeError):
        current_def.calculate_current(
            test_args.Q0, test_args.Q1, 100)

