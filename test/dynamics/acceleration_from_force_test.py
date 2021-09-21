from symplyphysics import (
    units, convert_to, SI
)
from symplyphysics.laws.dynamics import acceleration_from_force as newton_law2

def test_basic_force():
    m = units.Quantity('m')
    SI.set_quantity_dimension(m, units.mass)
    SI.set_quantity_scale_factor(m, 1 * units.kilogram)
    a = units.Quantity('a')
    SI.set_quantity_dimension(a, units.acceleration)
    SI.set_quantity_scale_factor(a, 3 * units.meter / units.second**2)

    result = newton_law2.calculate_force(m, a)
    result_force = convert_to(result, units.newton).subs(units.newton, 1).evalf(2)
    expected_mass = convert_to(m, units.kilogram).subs(units.kilogram, 1).evalf(2)
    expected_acceleration = convert_to(a, units.meter / units.second**2).subs({
            units.meter: 1, units.second: 1}).evalf(2)

    assert result_force == 3.0
    assert expected_mass == 1
    assert expected_acceleration == 3