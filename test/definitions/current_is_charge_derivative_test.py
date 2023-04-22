from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import current_is_charge_derivative as current_def


@fixture
def test_args():
    Q0 = Quantity(0 * units.coulomb)
    Q1 = Quantity(20 * units.coulomb)
    I0 = Quantity(0.5 * units.ampere)
    I1 = Quantity(0.6 * units.ampere)
    t = Quantity(5 * units.second)
    Args = namedtuple("Args", ["Q0", "Q1", "I0", "I1", "t"])
    return Args(Q0=Q0, Q1=Q1, I0=I0, I1=I1, t=t)


def test_basic_current(test_args):
    result = current_def.calculate_current(test_args.Q0, test_args.Q1, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.current)
    result_current = convert_to(result, current_def.definition_units_SI).subs(units.ampere,
        1).evalf(2)
    assert result_current == approx(4, 0.01)


def test_current_with_bad_charge(test_args):
    Qb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        current_def.calculate_current(Qb, test_args.Q1, test_args.t)
    with raises(errors.UnitsError):
        current_def.calculate_current(test_args.Q0, Qb, test_args.t)
    with raises(TypeError):
        current_def.calculate_current(100, test_args.Q1, test_args.t)
    with raises(TypeError):
        current_def.calculate_current(test_args.Q0, 100, test_args.t)


def test_current_with_bad_time(test_args):
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        current_def.calculate_current(test_args.Q0, test_args.Q1, tb)
    with raises(TypeError):
        current_def.calculate_current(test_args.Q0, test_args.Q1, 100)
