from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.euler_relations import gibbs_energy_formula

# Description
## The Gibbs energy of a one-component system consisting of 1e20 particles and having a chemical potential
## mu = -2 eV, is G = -32 J.

Args = namedtuple("Args", "mu n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mu = Quantity(-2 * units.electronvolt)
    n = 1e20
    return Args(mu=mu, n=n)


def test_law(test_args: Args) -> None:
    result = gibbs_energy_formula.calculate_gibbs_energy(test_args.mu, test_args.n)
    assert_equal(result, -32 * units.joule, tolerance=2e-3)


def test_bad_chemical_potential(test_args: Args) -> None:
    mub = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_energy_formula.calculate_gibbs_energy(mub, test_args.n)
    with raises(TypeError):
        gibbs_energy_formula.calculate_gibbs_energy(100, test_args.n)


def test_bad_particle_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_energy_formula.calculate_gibbs_energy(test_args.mu, nb)
