from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.thermodynamics import heat_from_sum_inner_energy_and_system_work as first_law_of_term


@fixture(name="test_args")
def fixture_test_args():
    W = Quantity(50 * units.joule)
    dU = Quantity(100 * units.joule)
    Args = namedtuple("Args", ["dU", "W"])
    return Args(dU=dU, W=W)


def test_basic_heat(test_args):
    result = first_law_of_term.calculate_heat(test_args.W, test_args.dU)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_heat = convert_to(result, units.joule).evalf(2)
    assert result_heat == approx(150, 0.01)


def test_bad_delta_inner_energy(test_args):
    error_dU = Quantity(test_args.dU * units.watt)
    with raises(errors.UnitsError):
        first_law_of_term.calculate_heat(error_dU, test_args.W)
    with raises(TypeError):
        first_law_of_term.calculate_heat(100, test_args.W)


def test_bad_work(test_args):
    error_W = Quantity(test_args.W * units.watt)
    with raises(errors.UnitsError):
        first_law_of_term.calculate_heat(test_args.dU, error_W)
    with raises(TypeError):
        first_law_of_term.calculate_heat(test_args.dU, 100)
