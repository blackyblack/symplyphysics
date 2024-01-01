from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.thermodynamics import first_law_of_termodynamics as first_law_of_term


@fixture(name='test_args')
def fixture_test_args():
    W = Quantity(50 * 10 ** 3 * units.joule)
    Q = Quantity(100 * 10 ** 3 * units.joule)
    Args = namedtuple('Args', ['Q', 'W'])
    return Args(Q=Q, W=W)


def test_basic_result(test_args):
    result = first_law_of_term.calculate_inner_energy(test_args.Q, test_args.W)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_delta_inner_energy = convert_to(result, units.joule).evalf(2)
    assert result_delta_inner_energy == approx(150 * 10 ** 3, 0.01)


def test_error_heat(test_args):
    error_q = Quantity(1 * units.watt)
    with raises(errors.UnitsError):
        first_law_of_term.calculate_inner_energy(error_q, test_args.W)
    with raises(TypeError):
        first_law_of_term.calculate_inner_energy(100, test_args.W)


def test_error_work(test_args):
    error_w = Quantity(1 * units.watt)
    with raises(errors.UnitsError):
        first_law_of_term.calculate_inner_energy(test_args.Q, error_w)
    with raises(TypeError):
        first_law_of_term.calculate_inner_energy(test_args.Q, 100)