from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic.rotational_inertia import rotational_inertia_of_particle as rotational_inertia_def

# Description
## A particle of mass m = 1 g rotates around an axis at a radius r = 3 m. Its rotational inertia
## should amount to 0.009 kg*m^2.

Args = namedtuple("Args", "m r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(1.0 * units.gram)
    r = Quantity(3.0 * units.meter)
    return Args(m=m, r=r)


def test_basic_law(test_args: Args) -> None:
    result = rotational_inertia_def.calculate_rotational_inertia(test_args.m, test_args.r)
    assert_equal(result, 9e-3 * units.kilogram * units.meter**2)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        rotational_inertia_def.calculate_rotational_inertia(mb, test_args.r)
    with raises(TypeError):
        rotational_inertia_def.calculate_rotational_inertia(100, test_args.r)


def test_bad_radius(test_args: Args) -> None:
    rb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        rotational_inertia_def.calculate_rotational_inertia(test_args.m, rb)
    with raises(TypeError):
        rotational_inertia_def.calculate_rotational_inertia(test_args.m, 100)
