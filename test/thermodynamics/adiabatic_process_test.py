from pytest import approx, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.thermodynamics import adiabatic_process

def test_basic_adiabatic_coefficient():
    P0 = units.Quantity('P0')
    SI.set_quantity_dimension(P0, units.pressure)
    SI.set_quantity_scale_factor(P0, 1 * units.pascal)
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.volume)
    SI.set_quantity_scale_factor(V0, 10 * units.liter)
    # Choose specific heats ratio
    y = 1.66

    result = adiabatic_process.calculate_adiabatic_coefficient(P0, V0, y)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure * units.volume**y)

    result_coeff = convert_to(result, units.pascal * units.liter**y).subs(units.liter, 1).subs(units.pascal, 1).evalf(4)
    assert result_coeff == approx(45.71, 0.01)
