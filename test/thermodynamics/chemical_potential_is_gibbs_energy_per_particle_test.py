from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    chemical_potential_is_gibbs_energy_per_particle as chemical_potential_law,)

# Description
## In a system of 100 particles, the Gibbs energy is 1000 J. The chemical potential
## of the system is therefore 10 J.

Args = namedtuple("Args", "g n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    g = Quantity(1000 * units.joule)
    n = 100
    return Args(g=g, n=n)


def test_law(test_args: Args) -> None:
    result = chemical_potential_law.calculate_chemical_potential(test_args.g, test_args.n)
    assert_equal(result, 10 * units.joule)


def test_bad_energy(test_args: Args) -> None:
    gb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(gb, test_args.n)
    with raises(TypeError):
        chemical_potential_law.calculate_chemical_potential(100, test_args.n)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.g, nb)
