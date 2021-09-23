from pytest import approx, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration

def test_basic_acceleration():
    v0 = units.Quantity('v0')
    SI.set_quantity_dimension(v0, units.velocity)
    SI.set_quantity_scale_factor(v0, 1 * units.meter / units.second)
    v1 = units.Quantity('v1')
    SI.set_quantity_dimension(v1, units.velocity)
    SI.set_quantity_scale_factor(v1, 20 * units.meter / units.second)
    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.time)
    SI.set_quantity_scale_factor(t, 5 * units.second)

    result = acceleration.calculate_linear_acceleration(v0, v1, t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.acceleration)

    result_acceleration = convert_to(result, acceleration.definition_dimension_SI).subs({
            units.meter: 1, units.second: 1}).evalf(2)
    assert result_acceleration == approx(3.8, 0.01)

def test_bad_velocity():
    v0 = units.Quantity('v0')
    v1 = units.Quantity('v1')
    SI.set_quantity_dimension(v1, units.velocity)
    SI.set_quantity_scale_factor(v1, 20 * units.meter / units.second)
    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.time)
    SI.set_quantity_scale_factor(t, 5 * units.second)

    # Make v0 invalid
    SI.set_quantity_dimension(v0, units.length)
    SI.set_quantity_scale_factor(v0, 1 * units.meter)
    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(v0, v1, t)

    # Make v0 valid and v1 invalid
    SI.set_quantity_dimension(v0, units.velocity)
    SI.set_quantity_scale_factor(v0, 1 * units.meter / units.second)
    SI.set_quantity_dimension(v1, units.length)
    SI.set_quantity_scale_factor(v1, 1 * units.meter)
    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(v0, v1, t)

    # Make v0 and v1 valid
    SI.set_quantity_dimension(v1, units.velocity)
    SI.set_quantity_scale_factor(v1, 20 * units.meter / units.second)

    with raises(TypeError):
        acceleration.calculate_linear_acceleration(100, v1, t)

    with raises(TypeError):
        acceleration.calculate_linear_acceleration(v0, 100, t)

def test_bad_time():
    v0 = units.Quantity('v0')
    SI.set_quantity_dimension(v0, units.velocity)
    SI.set_quantity_scale_factor(v0, 1 * units.meter / units.second)
    v1 = units.Quantity('v1')
    SI.set_quantity_dimension(v1, units.velocity)
    SI.set_quantity_scale_factor(v1, 20 * units.meter / units.second)
    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.length)
    SI.set_quantity_scale_factor(t, 1 * units.meter)

    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(v0, v1, t)
    
    with raises(TypeError):
        acceleration.calculate_linear_acceleration(v0, v1, 100)
