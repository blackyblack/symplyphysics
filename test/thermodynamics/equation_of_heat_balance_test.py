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
from symplyphysics.laws.thermodynamics import equation_of_heat_balance as heat_balance_law


@fixture(name="test_args")
def test_args_fixture():
    q1 = Quantity(3 * units.joules)
    q2 = Quantity(-5 * units.joules)
    Args = namedtuple("Args", ["q1", "q2"])
    return Args(q1=q1, q2=q2)


def test_basic_current(test_args):
    result = heat_balance_law.calculate_amount_energy([test_args.q1])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_current = convert_to(result, units.joule).evalf(2)
    assert result_current == approx(-3, 0.01)


def test_three_current_array(test_args):
    result = heat_balance_law.calculate_amount_energy([test_args.q1, test_args.q2])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_current = convert_to(result, units.joule).evalf(2)
    assert result_current == approx(2, 0.01)


def test_array_empty():
    result = heat_balance_law.calculate_amount_energy([])
    assert SI.get_dimension_system().is_dimensionless(result.dimension)
    assert int(convert_to(result, S.One).n()) == 0


def test_array_bad_current(test_args):
    qb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        heat_balance_law.calculate_amount_energy([test_args.q1, qb])
    with raises(TypeError):
        heat_balance_law.calculate_amount_energy([test_args.q1, 100])
    with raises(errors.UnitsError):
        heat_balance_law.calculate_amount_energy([qb, test_args.q2])
    with raises(TypeError):
        heat_balance_law.calculate_amount_energy([100, test_args.q2])
    with raises(errors.UnitsError):
        heat_balance_law.calculate_amount_energy([qb, qb])
    with raises(TypeError):
        heat_balance_law.calculate_amount_energy([100, 100])
    with raises(TypeError):
        heat_balance_law.calculate_amount_energy(test_args.q1)
