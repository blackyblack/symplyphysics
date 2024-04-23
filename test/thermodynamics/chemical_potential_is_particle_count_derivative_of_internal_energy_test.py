from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    chemical_potential_is_particle_count_derivative_of_internal_energy as chemical_potential_law,)

# Description
## In a thermodynamic system, when the particle count changed from 100 to 102 particles, the internal
## energy lowered from 100 J to 80 J. The chemical potential of the system is therefore -10 J.

Args = namedtuple("Args", "n0 n1 u0 u1")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n0 = 100
    n1 = 102
    u0 = Quantity(100 * units.joule)
    u1 = Quantity(80 * units.joule)
    return Args(n0=n0, n1=n1, u0=u0, u1=u1)


def test_law(test_args: Args) -> None:
    result = chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1,
        test_args.u0, test_args.u1)
    assert_equal(result, -10 * units.joule)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(nb, test_args.n1, test_args.u0,
            test_args.u1)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, nb, test_args.u0,
            test_args.u1)


def test_bad_energy(test_args: Args) -> None:
    ub = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1, ub,
            test_args.u1)
    with raises(TypeError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1, 100,
            test_args.u1)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1,
            test_args.u0, ub)
    with raises(TypeError):
        chemical_potential_law.calculate_chemical_potential(test_args.n0, test_args.n1,
            test_args.u0, 100)
