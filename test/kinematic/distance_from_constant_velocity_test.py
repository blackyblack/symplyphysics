from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import distance_from_constant_velocity as movement_law

# Description
## If an object is traveling at a constant velocity of 12m/s, then calculate the distance covered by the object after 1 minute.


@fixture(name="test_args")
def test_args_fixture():
    x0 = Quantity(0 * units.meter)
    v = Quantity(12 * units.meter / units.second)
    t = Quantity(60 * units.second)
    Args = namedtuple("Args", ["x0", "v", "t"])
    return Args(x0=x0, v=v, t=t)


def test_basic_distance(test_args):
    result = movement_law.calculate_distance(test_args.x0, test_args.v, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_vector = convert_to(result, units.meter).evalf(2)
    assert result_vector == approx(720, 0.01)


def test_bad_distance(test_args):
    xb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(xb, test_args.v, test_args.t)
    with raises(TypeError):
        movement_law.calculate_distance(100, test_args.v, test_args.t)


def test_bad_velocity(test_args):
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(test_args.x0, vb, test_args.t)
    with raises(TypeError):
        movement_law.calculate_distance(test_args.x0, 100, test_args.t)


def test_bad_time(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        movement_law.calculate_distance(test_args.x0, test_args.v, tb)
    with raises(TypeError):
        movement_law.calculate_distance(test_args.x0, test_args.v, 100)
