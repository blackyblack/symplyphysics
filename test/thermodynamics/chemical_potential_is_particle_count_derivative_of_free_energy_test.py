from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    chemical_potential_is_particle_count_derivative_of_free_energy as chemical_potential_law,)

# Description
## In a thermodynamic system, when the particle count changed from 100 to 101 particles, the free energy
## dropped from 100 J to 80 J. The chemical potential of the system is therefore -20 J.

Args = namedtuple("Args", "n0 n1 f0 f1")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n0 = 100
    n1 = 101
    f0 = Quantity(100 * units.joule)
    f1 = Quantity(80 * units.joule)
    return Args(n0=n0, n1=n1, f0=f0, f1=f1)


def test_law(test_args: Args) -> None:
    result = chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1,
        test_args.f0, test_args.f1)
    assert_equal(result, -20 * units.joule)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(nb, test_args.n1, test_args.f0,
            test_args.f1)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, nb, test_args.f0,
            test_args.f1)


def test_bad_energy(test_args: Args) -> None:
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1, fb,
            test_args.f1)
    with raises(TypeError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1, 100,
            test_args.f1)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1,
            test_args.f0, fb)
    with raises(TypeError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1,
            test_args.f0, 100)
