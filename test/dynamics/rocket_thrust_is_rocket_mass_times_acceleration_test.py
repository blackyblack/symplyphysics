from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import rocket_thrust_is_rocket_mass_times_acceleration as rocket_law

# Description
## A rocket whose initial mass is 850 kg consumes fuel at the rate of 2.3 kg/s. Its acceleration
## is 7.6 m/s**2. Then, the rocket's speed relative to the exhaust gases is about 2.8 km/s.


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(850.0 * units.kilogram)
    r = Quantity(2.3 * units.kilogram / units.second)
    a = Quantity(7.6 * units.meter / units.second**2)
    Args = namedtuple("Args", "m r a")
    return Args(m=m, r=r, a=a)


def test_basic_law(test_args):
    result = rocket_law.calculate_relative_velocity(test_args.r, test_args.m, test_args.a)
    assert_equal(result, 2.809 * units.kilometer / units.second)


def test_bad_consumption_rate(test_args):
    rb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        rocket_law.calculate_relative_velocity(rb, test_args.m, test_args.a)
    with raises(TypeError):
        rocket_law.calculate_relative_velocity(100, test_args.m, test_args.a)


def test_bad_mass(test_args):
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        rocket_law.calculate_relative_velocity(test_args.r, mb, test_args.a)
    with raises(TypeError):
        rocket_law.calculate_relative_velocity(test_args.r, 100, test_args.a)


def test_bad_acceleration(test_args):
    ab = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        rocket_law.calculate_relative_velocity(test_args.r, test_args.m, ab)
    with raises(TypeError):
        rocket_law.calculate_relative_velocity(test_args.r, test_args.m, 100)
