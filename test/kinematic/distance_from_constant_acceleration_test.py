from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import distance_from_constant_acceleration as accelerated_distance_law

# Description
## Test example: first example from https://ege-study.ru/materialy-ege/kurs-fiziki-teoriya/zadachi-na-ravnouskorennoe-dvizhenie/
## Let the distance be counted from the point where braking began. Then S0 = 0 m.
## The initial braking speed is 15 m/s. Braking occurs with an acceleration of -0.2 m/(s^2), because the acceleration is directed against the motion.
## Then the braking time is equal 50 seconds and the braking distance will be equal 500 m


@fixture(name="test_args")
def test_args_fixture():
    s0 = Quantity(0 * units.meter)
    v0 = Quantity(15 * units.meter / units.second)
    a = Quantity(-0.2 * units.meter / units.second**2)
    t = Quantity(50 * units.second)
    Args = namedtuple("Args", ["s0", "v0", "a", "t"])
    return Args(s0=s0, v0=v0, a=a, t=t)


def test_basic_velocity(test_args):
    result = accelerated_distance_law.calculate_distance(test_args.s0, test_args.v0, test_args.a, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_velocity = convert_to(result, units.meters).evalf(5)
    assert result_velocity == approx(500, 0.01)


def test_bad_distance(test_args):
    sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        accelerated_distance_law.calculate_distance(sb, test_args.v0, test_args.a, test_args.t)
    with raises(TypeError):
        accelerated_distance_law.calculate_distance(100, test_args.v0, test_args.a, test_args.t)


def test_bad_velocity(test_args):
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        accelerated_distance_law.calculate_distance(test_args.s0, vb, test_args.a, test_args.t)
    with raises(TypeError):
        accelerated_distance_law.calculate_distance(test_args.s0, 100, test_args.a, test_args.t)


def test_bad_acceleration(test_args):
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        accelerated_distance_law.calculate_distance(test_args.s0, test_args.v0, ab, test_args.t)
    with raises(TypeError):
        accelerated_distance_law.calculate_distance(test_args.s0, test_args.v0, 100, test_args.t)


def test_bad_time(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        accelerated_distance_law.calculate_distance(test_args.s0, test_args.v0, test_args.a, tb)
    with raises(TypeError):
        accelerated_distance_law.calculate_distance(test_args.s0, test_args.v0, test_args.a, 100)
