from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors, units, convert_to, Quantity, SI,
)
from symplyphysics.definitions import velocity_is_movement_derivative as velocity_def

# Description
## Assume object starts moving with some velocity. After 5 seconds it's distance is 80 meters.
## Velocity should be 16 m/s.

@fixture
def test_args():
    S0 = Quantity(0 * units.meter)
    S1 = Quantity(80 * units.meters)    
    t = Quantity(5 * units.second)
    Args = namedtuple("Args", ["S0", "S1", "t"])
    return Args(S0=S0, S1=S1, t=t)

def test_basic_velocity(test_args):
    result = velocity_def.calculate_velocity(test_args.S0, test_args.S1, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity)
    result_current = convert_to(result, velocity_def.definition_units_SI).subs({
        units.meter: 1, units.second: 1}).evalf(2)
    assert result_current == approx(16, 0.01)

def test_velocity_with_bad_distance(test_args):
    Sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_def.calculate_velocity(Sb, test_args.S1, test_args.t)
    with raises(errors.UnitsError):
        velocity_def.calculate_velocity(test_args.S0, Sb, test_args.t)
    with raises(TypeError):
        velocity_def.calculate_velocity(100, test_args.S1, test_args.t)
    with raises(TypeError):
        velocity_def.calculate_velocity(test_args.S0, 100, test_args.t)
    
def test_velocity_with_bad_time(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_def.calculate_velocity(test_args.S0, test_args.S1, tb)
    with raises(TypeError):
        velocity_def.calculate_velocity(test_args.S0, test_args.S1, 100)
