from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as movement_law

# Description
## The man starts running with initial speed 10m/s and is getting tired while running (decreases his velocity) with -0.1 m/s/s. How far this man can run away in 10 seconds?


@fixture(name="test_args")
def test_args_fixture():
    V0 = Quantity(10 * units.meter / units.second)
    a = Quantity(-0.1 * units.meter / units.second**2)
    t = Quantity(10 * units.second)
    Args = namedtuple("Args", ["V0", "a", "t"])
    return Args(V0=V0, a=a, t=t)


def test_basic_distance(test_args):
    result = movement_law.calculate_distance(test_args.V0, test_args.a, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_vector = convert_to(result, units.meter).evalf(2)
    assert result_vector == approx(95, 0.01)


def test_bad_velocity(test_args):
    Vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(Vb, test_args.a, test_args.t)
    with raises(TypeError):
        movement_law.calculate_distance(100, test_args.a, test_args.t)


def test_bad_acceleration(test_args):
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(test_args.V0, ab, test_args.t)
    with raises(TypeError):
        movement_law.calculate_distance(test_args.V0, 100, test_args.t)


def test_bad_time(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(test_args.V0, test_args.a, tb)
    with raises(TypeError):
        movement_law.calculate_distance(test_args.V0, test_args.a, 100)
