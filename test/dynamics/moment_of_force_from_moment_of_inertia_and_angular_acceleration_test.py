from collections import namedtuple
from pytest import approx, fixture, raises
from sympy.physics.units.systems import SI

from symplyphysics import (
    errors,
    units,
    Quantity,
    convert_to,
)
from symplyphysics.laws.dynamics import moment_of_force_from_moment_of_inertia_and_angular_acceleration as moment_force_law


@fixture(name="test_args")
def test_args_fixture():
    epsilon = Quantity(1 * units.radians / (units.second**2))
    i = Quantity(3 * units.kilograms * units.meters**2)
    Args = namedtuple("Args", ["i", "epsilon"])
    return Args(i=i, epsilon=epsilon)


def test_basic_moment_of_force(test_args):
    result = moment_force_law.calculate_moment_of_force(test_args.i, test_args.epsilon)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force * units.length)
    result_moment_of_force = convert_to(result, units.newtons * units.meters).evalf(2)
    assert result_moment_of_force == approx(3, 0.01)


def test_bad_moment_of_inertia(test_args):
    ib = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        moment_force_law.calculate_moment_of_force(ib, test_args.epsilon)
    with raises(TypeError):
        moment_force_law.calculate_moment_of_force(100, test_args.i)


def test_bad_angular_acceleration(test_args):
    epsilon_b = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        moment_force_law.calculate_moment_of_force(test_args.i, epsilon_b)
    with raises(TypeError):
        moment_force_law.calculate_moment_of_force(test_args.i, 100)
