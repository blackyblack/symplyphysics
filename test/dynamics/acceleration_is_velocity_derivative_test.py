#!/usr/bin/env python3
from pytest import approx

from symplyphysics import (
    units, convert_to, SI
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
    result_acceleration = convert_to(result, acceleration.definition_dimension_SI).subs({
            units.meter: 1, units.second: 1}).evalf(2)
    expected_initial_speed = convert_to(v0, units.meter / units.second).subs({
            units.meter: 1, units.second: 1}).evalf(2)
    expected_final_speed = convert_to(v1, units.meter / units.second).subs({
            units.meter: 1, units.second: 1}).evalf(2)
    expected_time = convert_to(t, units.second).subs(units.second, 1).evalf(2)

    assert result_acceleration == approx(3.8, 0.01)
    assert expected_initial_speed == 1
    assert expected_final_speed == 20
    assert expected_time == 5
