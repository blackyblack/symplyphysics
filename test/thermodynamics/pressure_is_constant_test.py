from pytest import approx, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.thermodynamics import pressure_is_constant as gay_lussacs_law

def test_basic_volume():
    t0 = units.Quantity('t0')
    SI.set_quantity_dimension(t0, units.temperature)
    SI.set_quantity_scale_factor(t0, 1 * units.kelvin)
    t1 = units.Quantity('t1')
    SI.set_quantity_dimension(t1, units.temperature)
    SI.set_quantity_scale_factor(t1, 2 * units.kelvin)
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.volume)
    SI.set_quantity_scale_factor(V0, 1 * units.liter)

    result = gay_lussacs_law.calculate_volume(t0, V0, t1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.volume)

    result_volume = convert_to(result, units.liter).subs(units.liter, 1).evalf(2)
    assert result_volume == approx(2.0, 0.01)

def test_bad_temperature():
    t0 = units.Quantity('t0')
    t1 = units.Quantity('t1')
    SI.set_quantity_dimension(t1, units.temperature)
    SI.set_quantity_scale_factor(t1, 2 * units.kelvin)
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.volume)
    SI.set_quantity_scale_factor(V0, 1 * units.liter)

    # Make t0 invalid
    SI.set_quantity_dimension(t0, units.length)
    SI.set_quantity_scale_factor(t0, 1 * units.meter)
    with raises(errors.UnitsError):
        gay_lussacs_law.calculate_volume(t0, V0, t1)

    # Make t0 valid and t1 invalid
    SI.set_quantity_dimension(t0, units.temperature)
    SI.set_quantity_scale_factor(t0, 1 * units.kelvin)
    SI.set_quantity_dimension(t1, units.length)
    SI.set_quantity_scale_factor(t1, 1 * units.meter)
    with raises(errors.UnitsError):
        gay_lussacs_law.calculate_volume(t0, V0, t1)

    # Make t0 and t1 valid
    SI.set_quantity_dimension(t0, units.temperature)
    SI.set_quantity_scale_factor(t0, 1 * units.kelvin)
    SI.set_quantity_dimension(t1, units.temperature)
    SI.set_quantity_scale_factor(t1, 2 * units.kelvin)

    with raises(TypeError):
        gay_lussacs_law.calculate_volume(100, V0, t1)
    
    with raises(TypeError):
        gay_lussacs_law.calculate_volume(t0, V0, 100)

def test_bad_volume():
    t0 = units.Quantity('t0')
    SI.set_quantity_dimension(t0, units.temperature)
    SI.set_quantity_scale_factor(t0, 1 * units.kelvin)
    t1 = units.Quantity('t1')
    SI.set_quantity_dimension(t1, units.temperature)
    SI.set_quantity_scale_factor(t1, 2 * units.kelvin)
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.length)
    SI.set_quantity_scale_factor(V0, 1 * units.meter)

    with raises(errors.UnitsError):
        gay_lussacs_law.calculate_volume(t0, V0, t1)

    with raises(TypeError):
        gay_lussacs_law.calculate_volume(t0, 100, t1)
