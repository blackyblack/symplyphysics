from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.statistical_physics.classical_statistics.maxwell_boltzmann import maxwell_boltzmann_most_probable_speed as law

# Description
## The most probable speed of Argon particles (m = 39.948 u) at equilibrium temperature T = 150 K
## is 250 m/s.

Args = namedtuple("Args", "t m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(150 * units.kelvin)
    m = Quantity(39.948 * units.amu)
    return Args(t=t, m=m)


def test_law(test_args: Args) -> None:
    result = law.calculate_most_probable_speed(test_args.t, test_args.m)
    assert_equal(result, 250 * units.meter / units.second)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_most_probable_speed(tb, test_args.m)
    with raises(TypeError):
        law.calculate_most_probable_speed(100, test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_most_probable_speed(test_args.t, mb)
    with raises(TypeError):
        law.calculate_most_probable_speed(test_args.t, 100)
