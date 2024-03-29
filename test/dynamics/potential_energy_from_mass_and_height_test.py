from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import potential_energy_from_mass_and_height as potential_energy

# How much potential energy does the body with a mass of 9 grams have at
# an altitude of 500 meters?

Args = namedtuple("Args", ["m", "h"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(9 * units.gram)
    h = Quantity(500 * units.meter)
    return Args(m=m, h=h)


def test_basic_energy(test_args: Args) -> None:
    result = potential_energy.calculate_potential_energy(test_args.m, test_args.h)
    assert_equal(result, 44.13 * units.joule)


def test_bad_body_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_energy.calculate_potential_energy(mb, test_args.h)
    with raises(TypeError):
        potential_energy.calculate_potential_energy(100, test_args.h)


def test_bad_height(test_args: Args) -> None:
    hb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_energy.calculate_potential_energy(test_args.m, hb)
    with raises(TypeError):
        potential_energy.calculate_potential_energy(test_args.m, 100)
