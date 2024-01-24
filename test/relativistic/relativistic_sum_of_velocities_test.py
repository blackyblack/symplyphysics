from collections import namedtuple

from pytest import approx, fixture, raises
from sympy.physics.units import speed_of_light

from symplyphysics import Quantity, convert_to, errors, units
from symplyphysics.laws.relativistic import relativistic_sum_of_velocities

# Description
# v1 = 0.5 * c, v2 = 0.5 * c
# v = (v1 + v2) / (1 + (v1 * v2) / c**2) = (0.5 * c + 0.5 * c) / (1 + (0.25 * c**2) / c**2)
# = c / (1.25 * c**2) = 0.8 * c
# Expected result is 0.8 * c


@fixture(name="test_args")
def test_args_fixture():
    v1 = Quantity(0.5 * speed_of_light)
    v2 = Quantity(0.5 * speed_of_light)
    Args = namedtuple("Args", ["v1", "v2"])
    return Args(v1=v1, v2=v2)


def test_basic_sum(test_args):
    result = relativistic_sum_of_velocities.calculate_relativistic_sum_of_velocities(
        test_args.v1, test_args.v2)
    result_velocity = convert_to(result, units.meter / units.second).evalf(4)
    result_expected = convert_to(
        Quantity(0.8 * speed_of_light), units.meter / units.second).evalf(4)
    assert result_velocity == approx(result_expected, 0.01)


def test_bad_velocity(test_args):
    mv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relativistic_sum_of_velocities.calculate_relativistic_sum_of_velocities(
            mv, test_args.v2)
    with raises(TypeError):
        relativistic_sum_of_velocities.calculate_relativistic_sum_of_velocities(
            100, test_args.v1)
    with raises(errors.UnitsError):
        relativistic_sum_of_velocities.calculate_relativistic_sum_of_velocities(
            test_args.v1, mv)
    with raises(TypeError):
        relativistic_sum_of_velocities.calculate_relativistic_sum_of_velocities(
            test_args.v1, 100)
