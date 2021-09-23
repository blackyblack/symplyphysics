from pytest import approx, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.thermodynamics import pressure_from_temperature_and_volume as ideal_gas_law

# The volume of 1 mol of any gas at STP (Standard temperature, 273.15 K and pressure, 1 atm) is measured to be 22.414 L.
def test_basic_pressure():
    V = units.Quantity('V')
    SI.set_quantity_dimension(V, units.volume)
    SI.set_quantity_scale_factor(V, 22.414 * units.liter)
    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.temperature)
    SI.set_quantity_scale_factor(t, 273.15 * units.kelvin)
    n = units.Quantity('n')
    SI.set_quantity_dimension(n, units.amount_of_substance)
    SI.set_quantity_scale_factor(n, 1 * units.mole)

    result = ideal_gas_law.calculate_pressure(V, t, n)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)

    result_pressure = convert_to(result, units.pascal).subs(units.pascal, 1).evalf(2)
    assert result_pressure == approx(101.325e+3, 0.01)

    # Also check that calculated pressure = 1 atmosphere
    result_pressure_atms = convert_to(result, units.atm).subs(units.atm, 1).evalf(2)
    assert result_pressure_atms == approx(1.0, 0.01)

def test_bad_volume():
    V = units.Quantity('V')
    SI.set_quantity_dimension(V, units.length)
    SI.set_quantity_scale_factor(V, 1 * units.meter)
    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.temperature)
    SI.set_quantity_scale_factor(t, 273.15 * units.kelvin)
    n = units.Quantity('n')
    SI.set_quantity_dimension(n, units.amount_of_substance)
    SI.set_quantity_scale_factor(n, 1 * units.mole)

    with raises(errors.UnitsError):
        ideal_gas_law.calculate_pressure(V, t, n)

    with raises(TypeError):
        ideal_gas_law.calculate_pressure(100, t, n)

def test_bad_temperature():
    V = units.Quantity('V')
    SI.set_quantity_dimension(V, units.volume)
    SI.set_quantity_scale_factor(V, 22.414 * units.liter)
    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.length)
    SI.set_quantity_scale_factor(t, 1 * units.meter)
    n = units.Quantity('n')
    SI.set_quantity_dimension(n, units.amount_of_substance)
    SI.set_quantity_scale_factor(n, 1 * units.mole)

    with raises(errors.UnitsError):
        ideal_gas_law.calculate_pressure(V, t, n)

    with raises(TypeError):
        ideal_gas_law.calculate_pressure(V, 100, n)

def test_bad_mole_count():
    V = units.Quantity('V')
    SI.set_quantity_dimension(V, units.volume)
    SI.set_quantity_scale_factor(V, 22.414 * units.liter)
    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.temperature)
    SI.set_quantity_scale_factor(t, 273.15 * units.kelvin)
    n = units.Quantity('n')
    SI.set_quantity_dimension(n, units.length)
    SI.set_quantity_scale_factor(n, 1 * units.meter)

    with raises(errors.UnitsError):
        ideal_gas_law.calculate_pressure(V, t, n)

    with raises(TypeError):
        ideal_gas_law.calculate_pressure(V, t, 100)
