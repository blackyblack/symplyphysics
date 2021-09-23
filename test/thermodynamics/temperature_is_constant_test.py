from pytest import approx, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.thermodynamics import temperature_is_constant as boyles_law

def test_basic_volume():
    P0 = units.Quantity('P0')
    SI.set_quantity_dimension(P0, units.pressure)
    SI.set_quantity_scale_factor(P0, 1 * units.pascal)
    P1 = units.Quantity('P1')
    SI.set_quantity_dimension(P1, units.pressure)
    SI.set_quantity_scale_factor(P1, 2 * units.pascal)
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.volume)
    SI.set_quantity_scale_factor(V0, 1 * units.liter)

    result = boyles_law.calculate_volume(P0, V0, P1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.volume)

    result_volume = convert_to(result, units.liter).subs(units.liter, 1).evalf(2)
    assert result_volume == approx(0.5, 0.01)

def test_bad_pressure():
    P0 = units.Quantity('P0')
    P1 = units.Quantity('P1')
    SI.set_quantity_dimension(P1, units.pressure)
    SI.set_quantity_scale_factor(P1, 2 * units.pascal)
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.volume)
    SI.set_quantity_scale_factor(V0, 1 * units.liter)

    # Make P0 invalid
    SI.set_quantity_dimension(P0, units.length)
    SI.set_quantity_scale_factor(P0, 1 * units.meter)
    with raises(errors.UnitsError):
        boyles_law.calculate_volume(P0, V0, P1)

    # Make P0 valid and P1 invalid
    SI.set_quantity_dimension(P0, units.pressure)
    SI.set_quantity_scale_factor(P0, 1 * units.pascal)
    SI.set_quantity_dimension(P1, units.length)
    SI.set_quantity_scale_factor(P1, 1 * units.meter)
    with raises(errors.UnitsError):
        boyles_law.calculate_volume(P0, V0, P1)

    # Make P0 and P1 valid
    SI.set_quantity_dimension(P0, units.pressure)
    SI.set_quantity_scale_factor(P0, 1 * units.pascal)
    SI.set_quantity_dimension(P1, units.pressure)
    SI.set_quantity_scale_factor(P1, 2 * units.pascal)

    with raises(TypeError):
        boyles_law.calculate_volume(100, V0, P1)
    
    with raises(TypeError):
        boyles_law.calculate_volume(P0, V0, 100)

def test_bad_volume():
    P0 = units.Quantity('P0')
    SI.set_quantity_dimension(P0, units.pressure)
    SI.set_quantity_scale_factor(P0, 1 * units.pascal)
    P1 = units.Quantity('P1')
    SI.set_quantity_dimension(P1, units.pressure)
    SI.set_quantity_scale_factor(P1, 2 * units.pascal)
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.length)
    SI.set_quantity_scale_factor(V0, 1 * units.meter)

    with raises(errors.UnitsError):
        boyles_law.calculate_volume(P0, V0, P1)

    with raises(TypeError):
        boyles_law.calculate_volume(P0, 100, P1)
