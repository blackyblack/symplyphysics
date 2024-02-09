from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity import gravitational_potential_energy

# Description
## Two particles, one with mass m1 = 4 g and the other with mass m2 = 10 g, are 10 cm apart
## from one another. The gravitational potential energy of this system is -2.67e-14 J.


@fixture(name="test_args")
def test_args_fixture():
    m1 = Quantity(4.0 * units.gram)
    m2 = Quantity(10.0 * units.gram)
    r = Quantity(10.0 * units.centimeter)
    Args = namedtuple("Args", "m1 m2 r")
    return Args(m1=m1, m2=m2, r=r)


def test_law(test_args):
    result = gravitational_potential_energy.calculate_gravitational_potential_energy(
        test_args.m1, test_args.m2, test_args.r)
    assert_equal(result, -2.67e-14 * units.joule)


def test_bad_masses(test_args):
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        gravitational_potential_energy.calculate_gravitational_potential_energy(
            mb, test_args.m2, test_args.r)
    with raises(TypeError):
        gravitational_potential_energy.calculate_gravitational_potential_energy(
            100, test_args.m2, test_args.r)
    with raises(errors.UnitsError):
        gravitational_potential_energy.calculate_gravitational_potential_energy(
            test_args.m1, mb, test_args.r)
    with raises(TypeError):
        gravitational_potential_energy.calculate_gravitational_potential_energy(
            test_args.m1, 100, test_args.r)


def test_bad_distance(test_args):
    rb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        gravitational_potential_energy.calculate_gravitational_potential_energy(
            test_args.m1, test_args.m2, rb)
    with raises(TypeError):
        gravitational_potential_energy.calculate_gravitational_potential_energy(
            test_args.m1, test_args.m2, 100)
