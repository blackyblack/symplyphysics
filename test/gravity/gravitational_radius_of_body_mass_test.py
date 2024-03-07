from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.gravity import gravitational_radius_of_body_mass as radius_law

# Description
## The mass of the Earth is 5.98e24 kilograms. Then its gravitational radius is 8.882e-3 meters.
## http://nuclphys.sinp.msu.ru/problems/095.htm

Args = namedtuple("Args", ["body_mass"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    body_mass = Quantity(5.98e24 * units.kilogram)

    return Args(body_mass=body_mass)


def test_basic_gravitational_radius(test_args: Args) -> None:
    result = radius_law.calculate_radius(test_args.body_mass)
    assert_equal(result, 8.882e-3 * units.meter)


def test_bad_body_mass(test_args: Args) -> None:
    body_mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(body_mass)
    with raises(TypeError):
        radius_law.calculate_radius(100)
