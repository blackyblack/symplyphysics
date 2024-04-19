from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics.euler_relations import enthalpy_formula

# Description
## The enthalpy of a system at temperature T = 100 K, having entropy S = 1 kJ/K and chemical potential
## mu = -1 eV and N = 1e22 particles is H = 98 kJ.

Args = namedtuple("Args", "t s mu n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(100 * units.kelvin)
    s = Quantity(1 * prefixes.kilo * units.joule / units.kelvin)
    mu = Quantity(-1 * units.electronvolt)
    n = 1e22
    return Args(t=t, s=s, mu=mu, n=n)


def test_law(test_args: Args) -> None:
    result = enthalpy_formula.calculate_enthalpy(test_args.t, test_args.s, test_args.mu,
        test_args.n)
    assert_equal(result, 98 * prefixes.kilo * units.joule, tolerance=5e-3)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_formula.calculate_enthalpy(tb, test_args.s, test_args.mu, test_args.n)
    with raises(TypeError):
        enthalpy_formula.calculate_enthalpy(100, test_args.s, test_args.mu, test_args.n)


def test_bad_entropy(test_args: Args) -> None:
    sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_formula.calculate_enthalpy(test_args.t, sb, test_args.mu, test_args.n)
    with raises(TypeError):
        enthalpy_formula.calculate_enthalpy(test_args.t, 100, test_args.mu, test_args.n)


def test_bad_chemical_potential(test_args: Args) -> None:
    mub = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_formula.calculate_enthalpy(test_args.t, test_args.s, mub, test_args.n)
    with raises(TypeError):
        enthalpy_formula.calculate_enthalpy(test_args.t, test_args.s, 100, test_args.n)


def test_bad_particle_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_formula.calculate_enthalpy(test_args.t, test_args.s, test_args.mu, nb)
