from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    quantities,
)
from symplyphysics.laws.thermodynamics import chemical_potential_of_ideal_gas as chemical_potential

# Description
## The chemical potential of Argon gas at standard conditions (n = 2.69e22 1/L, T = 273.15 K,
## lambda = 0.16 Å) is mu = -0.377 eV.

Args = namedtuple("Args", "t n l")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = quantities.standard_conditions_temperature
    n = Quantity(2.69e22 / units.liter)
    l = Quantity(0.16 * units.angstrom)
    return Args(t=t, n=n, l=l)


def test_law(test_args: Args) -> None:
    result = chemical_potential.calculate_chemical_potential(test_args.t, test_args.n, test_args.l)
    assert_equal(result, -0.377 * units.electronvolt)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential.calculate_chemical_potential(tb, test_args.n, test_args.l)
    with raises(TypeError):
        chemical_potential.calculate_chemical_potential(100, test_args.n, test_args.l)


def test_bad_concentration(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential.calculate_chemical_potential(test_args.t, nb, test_args.l)
    with raises(TypeError):
        chemical_potential.calculate_chemical_potential(test_args.t, 100, test_args.l)


def test_bad_wavelength(test_args: Args) -> None:
    lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        chemical_potential.calculate_chemical_potential(test_args.t, test_args.n, lb)
    with raises(TypeError):
        chemical_potential.calculate_chemical_potential(test_args.t, test_args.n, 100)
