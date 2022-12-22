from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors, pi
)
from symplyphysics.definitions import angular_velocity_is_angle_derivative as angular_velocity_def
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

# Description
## Assume object starts rotating with some angular velocity. After 5 seconds it rotates to 180 degrees.
## Velocity should be pi/5 radian/sec.

@fixture
def test_args():
    a0 = units.Quantity('a0')
    SI.set_quantity_dimension(a0, angle_type)
    SI.set_quantity_scale_factor(a0, 0 * units.radian)
    a1 = units.Quantity('a1')
    SI.set_quantity_dimension(a1, angle_type)
    SI.set_quantity_scale_factor(a1, pi * units.radian)    
    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.time)
    SI.set_quantity_scale_factor(t, 5 * units.second)

    Args = namedtuple('Args', ['a0', 'a1', 't'])
    return Args(a0=a0, a1=a1, t=t)


def test_basic_velocity(test_args):
    result = angular_velocity_def.calculate_velocity(test_args.a0, test_args.a1, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.time)
    result_velocity = convert_to(result, angular_velocity_def.definition_dimension_SI).subs({units.radian: 1, units.second: 1}).evalf(2)
    assert result_velocity == approx(0.63, 0.01)


def test_velocity_with_bad_angle(test_args):
    a0b = units.Quantity('a0b')
    SI.set_quantity_dimension(a0b, units.charge)
    SI.set_quantity_scale_factor(a0b, 1 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_def.calculate_velocity(
            a0b, test_args.a1, test_args.t)
    
    with raises(errors.UnitsError):
        angular_velocity_def.calculate_velocity(
            test_args.a0, a0b, test_args.t)           

    with raises(TypeError):
        angular_velocity_def.calculate_velocity(
            100, test_args.a1, test_args.t)

    with raises(TypeError):
        angular_velocity_def.calculate_velocity(
            test_args.a0, 100, test_args.t)

    
def test_velocity_with_bad_time(test_args):
    tb = units.Quantity('tb')
    SI.set_quantity_dimension(tb, units.charge)
    SI.set_quantity_scale_factor(tb, 1 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_def.calculate_velocity(
            test_args.a0, test_args.a1, tb)
    
    with raises(TypeError):
        angular_velocity_def.calculate_velocity(
            test_args.a0, test_args.a1, 100)

