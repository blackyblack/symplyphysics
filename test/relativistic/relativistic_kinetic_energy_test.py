from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.relativistic import relativistic_kinetic_energy as law

Args = namedtuple("Args", "g m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    g = 10
    m = units.electron_rest_mass
    return Args(g=g, m=m)


def test_law(test_args: Args) -> None:
    result = law.calculate_kinetic_energy(test_args.g, test_args.m)
    assert_equal(result, 4.6 * prefixes.mega * units.electronvolt)


def test_bad_factor(test_args: Args) -> None:
    gb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_kinetic_energy(gb, test_args.m)

    with raises(ValueError):
        law.calculate_kinetic_energy(0.3, test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_kinetic_energy(test_args.g, mb)
    with raises(TypeError):
        law.calculate_kinetic_energy(test_args.g, 100)
