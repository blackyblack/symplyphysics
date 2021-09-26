from pytest import approx, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.thermodynamics import zero_heat_transfer

def test_basic_pressure():
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

    result = zero_heat_transfer.calculate_pressure(n, t0, V0, V1, y)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)

    result_coeff = convert_to(result, units.pascal).subs(units.pascal, 1).evalf(8)
    assert result_coeff == approx(2631.02, 0.01)

def test_bad_mole_count():
    n = units.Quantity('n')
    SI.set_quantity_dimension(n, units.length)
    SI.set_quantity_scale_factor(n, 1 * units.meter)
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

    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(n, t0, V0, V1, y)

    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(100, t0, V0, V1, y)
