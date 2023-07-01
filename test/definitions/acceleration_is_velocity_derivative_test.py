from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration


@fixture(name="test_args")
def test_args_fixture():
    v0 = Quantity(1 * units.meter / units.second)
    v1 = Quantity(20 * units.meter / units.second)
    t = Quantity(5 * units.second)
    Args = namedtuple("Args", ["v0", "v1", "t"])
    return Args(v0=v0, v1=v1, t=t)


def test_basic_acceleration(test_args):
    result = acceleration.calculate_linear_acceleration(test_args.v0, test_args.v1, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.acceleration)
    result_acceleration = convert_to(result, acceleration.definition_units_SI).subs({
        units.meter: 1,
        units.second: 1
    }).evalf(2)
    assert result_acceleration == approx(3.8, 0.01)


def test_bad_velocity(test_args):
    vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(vb, test_args.v1, test_args.t)
    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(test_args.v0, vb, test_args.t)
    with raises(TypeError):
        acceleration.calculate_linear_acceleration(100, test_args.v1, test_args.t)
    with raises(TypeError):
        acceleration.calculate_linear_acceleration(test_args.v0, 100, test_args.t)


def test_bad_time(test_args):
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        acceleration.calculate_linear_acceleration(test_args.v0, test_args.v1, tb)
    with raises(TypeError):
        acceleration.calculate_linear_acceleration(test_args.v0, test_args.v1, 100)
