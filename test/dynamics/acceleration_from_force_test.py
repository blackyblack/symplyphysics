from pytest import approx, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.dynamics import acceleration_from_force as newton_second_law

def test_basic_force():
    m = units.Quantity('m')
    SI.set_quantity_dimension(m, units.mass)
    SI.set_quantity_scale_factor(m, 1 * units.kilogram)
    a = units.Quantity('a')
    SI.set_quantity_dimension(a, units.acceleration)
    SI.set_quantity_scale_factor(a, 3 * units.meter / units.second**2)

    result = newton_second_law.calculate_force(m, a)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)

    result_force = convert_to(result, units.newton).subs(units.newton, 1).evalf(2)
    assert result_force == approx(3.0, 0.01)

def test_bad_mass():
    m = units.Quantity('m')
    SI.set_quantity_dimension(m, units.length)
    SI.set_quantity_scale_factor(m, 1 * units.meter)
    a = units.Quantity('a')
    SI.set_quantity_dimension(a, units.acceleration)
    SI.set_quantity_scale_factor(a, 3 * units.meter / units.second**2)

    with raises(errors.UnitsError):
        newton_second_law.calculate_force(m, a)

    with raises(TypeError):
        newton_second_law.calculate_force(100, a)

def test_bad_acceleration():
    m = units.Quantity('m')
    SI.set_quantity_dimension(m, units.mass)
    SI.set_quantity_scale_factor(m, 1 * units.kilogram)
    a = units.Quantity('a')
    SI.set_quantity_dimension(a, units.length)
    SI.set_quantity_scale_factor(a, 3 * units.meter)

    with raises(errors.UnitsError):
        newton_second_law.calculate_force(m, a)

    with raises(TypeError):
        newton_second_law.calculate_force(m, 100)
