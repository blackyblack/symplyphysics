from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import period_of_physical_pendulum as pendulum_law

# Description
## The pendulum's rotational inertia is 0.05 kg*m**2, its mass is 150 g and the distance
## from the pendulum's center of mass to pivot point is 15 cm. The oscillation period
## of the pendulum is 2.99 s.


@fixture(name="test_args")
def test_args_fixture():
    I = Quantity(0.05 * units.kilogram * units.meter**2)
    m = Quantity(150.0 * units.gram)
    h = Quantity(15.0 * units.centimeter)
    Args = namedtuple("Args", "I m h")
    return Args(I=I, m=m, h=h)


def test_law(test_args):
    result = pendulum_law.calculate_period(test_args.I, test_args.m, test_args.h)
    assert_equal(result, 2.99 * units.second)


def test_bad_rotational_inertia(test_args):
    Ib = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        pendulum_law.calculate_period(Ib, test_args.m, test_args.h)
    with raises(TypeError):
        pendulum_law.calculate_period(100, test_args.m, test_args.h)


def test_bad_mass(test_args):
    mb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        pendulum_law.calculate_period(test_args.I, mb, test_args.h)
    with raises(TypeError):
        pendulum_law.calculate_period(test_args.I, 100, test_args.h)


def test_bad_distance(test_args):
    hb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        pendulum_law.calculate_period(test_args.I, test_args.m, hb)
    with raises(TypeError):
        pendulum_law.calculate_period(test_args.I, test_args.m, 100)
