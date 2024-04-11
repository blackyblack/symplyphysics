from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    chemical_potential_via_derivative_with_respect_to_particle_count as chemical_potential_law,
)

# Description
## A thermodynamic system is described by some thermodynamic potential whose initial value was X = 100 J
## at particle count N = 10, and final value X = 140 J at particle count N = 12. The chemical potential
## of the system is mu = 20 J.

Args = namedtuple("Args", "x0 x1 n0 n1")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    x0 = Quantity(100 * units.joule)
    x1 = Quantity(140 * units.joule)
    n0 = 10
    n1 = 12
    return Args(x0=x0, x1=x1, n0=n0, n1=n1)


def test_law(test_args: Args) -> None:
    result = chemical_potential_law.calculate_chemical_potential(test_args.x0, test_args.x1, test_args.n0, test_args.n1)
    assert_equal(result, 20 * units.joule)


def test_bad_energy(test_args: Args) -> None:
    xb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(xb, test_args.x1, test_args.n0, test_args.n1)
    with raises(TypeError):
        chemical_potential_law.calculate_chemical_potential(100, test_args.x1, test_args.n0, test_args.n1)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.x0, xb, test_args.n0, test_args.n1)
    with raises(TypeError):
        chemical_potential_law.calculate_chemical_potential(test_args.x0, 100, test_args.n0, test_args.n1)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.x0, test_args.x1, nb, test_args.n1)
    with raises(errors.UnitsError):
        chemical_potential_law.calculate_chemical_potential(test_args.x0, test_args.x1, test_args.n0, nb)
