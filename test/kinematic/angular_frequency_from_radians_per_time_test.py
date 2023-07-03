from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    angle_type,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import angular_frequency_from_radians_per_time as frequency_def

# Description
## Circle is rotating to 9.42 radians in 5.4 seconds. It should have angular frequency of 1.744 rad/s.


@fixture(name="test_args")
def test_args_fixture():
    rotation = Quantity(9.42 * units.radian, dimension=angle_type)
    t = Quantity(5.4 * units.second)
    Args = namedtuple("Args", ["N", "t"])
    return Args(N=rotation, t=t)


def test_basic_frequency(test_args):
    result = frequency_def.calculate_frequency(test_args.N, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.time)
    result_frequency = convert_to(result, units.radian / units.second).evalf(2)
    assert result_frequency == approx(1.744, 0.01)


def test_bad_time(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_def.calculate_frequency(test_args.N, tb)
    with raises(TypeError):
        frequency_def.calculate_frequency(test_args.N, 100)


def test_bad_radians(test_args):
    Nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_def.calculate_frequency(Nb, test_args.t)
