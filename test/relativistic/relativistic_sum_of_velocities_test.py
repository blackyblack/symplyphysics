from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_approx, Quantity, convert_to, errors, units
from symplyphysics.laws.relativistic import relativistic_sum_of_velocities

# Compared with example at https://openstax.org/books/college-physics-2e/pages/28-4-relativistic-addition-of-velocities


@fixture(name="test_args")
def test_args_fixture():
    v1 = Quantity(0.5 * units.speed_of_light)
    v2 = Quantity(1 * units.speed_of_light)
    Args = namedtuple("Args", ["v1", "v2"])
    return Args(v1=v1, v2=v2)


def test_basic_sum(test_args):
    result = relativistic_sum_of_velocities.calculate_velocity(test_args.v1, test_args.v2)
    result_velocity = convert_to(result, units.meter / units.second).evalf(4)
    result_expected = convert_to(Quantity(units.speed_of_light),
        units.meter / units.second).evalf(4)
    assert_approx(result_velocity, result_expected)


def test_bad_velocity(test_args):
    mv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relativistic_sum_of_velocities.calculate_velocity(mv, test_args.v2)
    with raises(TypeError):
        relativistic_sum_of_velocities.calculate_velocity(100, test_args.v1)
    with raises(errors.UnitsError):
        relativistic_sum_of_velocities.calculate_velocity(test_args.v1, mv)
    with raises(TypeError):
        relativistic_sum_of_velocities.calculate_velocity(test_args.v1, 100)
