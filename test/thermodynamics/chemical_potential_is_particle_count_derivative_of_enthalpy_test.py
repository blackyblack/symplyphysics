from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    chemical_potential_is_particle_count_derivative_of_enthalpy as chemical_potential_law,
)

# Description
## In a thermodynamic system, when the particle count changed from 100 to 101 particles, the enthalpy
## increased from 80 J to 100 J. The chemical potential of the system is therefore 20 J.

Args = namedtuple("Args", "n0 n1 h0 h1")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n0 = 100
    n1 = 101
    h0 = Quantity(80 * units.joule)
    h1 = Quantity(100 * units.joule)
    return Args(n0=n0, n1=n1, h0=h0, h1=h1)


def test_law(test_args: Args) -> None:
    result = chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1, test_args.h0, test_args.h1)
    assert_equal(result, 20 * units.joule)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(nb, test_args.n1, test_args.h0, test_args.h1)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, nb, test_args.h0, test_args.h1)


def test_bad_energy(test_args: Args) -> None:
    hb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1, hb, test_args.h1)
    with raises(TypeError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1, 100, test_args.h1)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1, test_args.h0, hb)
    with raises(TypeError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1, test_args.h0, 100)
