from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.thermodynamics import zero_heat_transfer

@fixture
def test_args():
    n = units.Quantity('n')
    SI.set_quantity_dimension(n, units.amount_of_substance)
    SI.set_quantity_scale_factor(n, 1 * units.mole)
    t0 = units.Quantity('t0')
    SI.set_quantity_dimension(t0, units.temperature)
    SI.set_quantity_scale_factor(t0, 1 * units.kelvin)
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.volume)
    SI.set_quantity_scale_factor(V0, 1 * units.liter)
    V1 = units.Quantity('V1')
    SI.set_quantity_dimension(V1, units.volume)
    SI.set_quantity_scale_factor(V1, 2 * units.liter)
    # Choose specific heats ratio
    y = 1.66

    Args = namedtuple('Args', ['n', 't0', 'V0', 'V1', 'y'])
    return Args(n = n, t0 = t0, V0 = V0, V1 = V1, y = y)

def test_basic_pressure(test_args):
    result = zero_heat_transfer.calculate_pressure(
        test_args.n, test_args.t0, test_args.V0, test_args.V1, test_args.y)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)

    result_coeff = convert_to(result, units.pascal).subs(units.pascal, 1).evalf(8)
    assert result_coeff == approx(2631.02, 0.01)

def test_bad_mole_count(test_args):
    nb = units.Quantity('nb')
    SI.set_quantity_dimension(nb, units.length)
    SI.set_quantity_scale_factor(nb, 1 * units.meter)

    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(nb, test_args.t0, test_args.V0, test_args.V1, test_args.y)

    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(100, test_args.t0, test_args.V0, test_args.V1, test_args.y)
