from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics.euler_relations import internal_energy_formula as formula

# Description
## The internal energy of a system at temperature T = 300 K, pressure p = 10 atm, with entropy S = 10 J/K,
## chemical potential mu = 100 J/K, and 1e6 particles, occupying 1 m**3 of space is U = 99 MJ.

Args = namedtuple("Args", "t s p v mu n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(300 * units.kelvin)
    s = Quantity(10 * units.joule / units.kelvin)
    p = Quantity(10 * units.atmosphere)
    v = Quantity(1 * units.meter**3)
    mu = Quantity(100 * units.joule)
    n = 1e6
    return Args(t=t, s=s, p=p, v=v, mu=mu, n=n)


def test_law(test_args: Args) -> None:
    result = formula.calculate_internal_energy(test_args.t, test_args.s, test_args.p, test_args.v, test_args.mu, test_args.n)
    assert_equal(result, 99 * prefixes.mega * units.joule)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        formula.calculate_internal_energy(tb, test_args.s, test_args.p, test_args.v, test_args.mu, test_args.n)
    with raises(TypeError):
        formula.calculate_internal_energy(100, test_args.s, test_args.p, test_args.v, test_args.mu, test_args.n)


def test_bad_entropy(test_args: Args) -> None:
    sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        formula.calculate_internal_energy(test_args.t, sb, test_args.p, test_args.v, test_args.mu, test_args.n)
    with raises(TypeError):
        formula.calculate_internal_energy(test_args.t, 100, test_args.p, test_args.v, test_args.mu, test_args.n)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        formula.calculate_internal_energy(test_args.t, test_args.s, pb, test_args.v, test_args.mu, test_args.n)
    with raises(TypeError):
        formula.calculate_internal_energy(test_args.t, test_args.s, 100, test_args.v, test_args.mu, test_args.n)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        formula.calculate_internal_energy(test_args.t, test_args.s, test_args.p, vb, test_args.mu, test_args.n)
    with raises(TypeError):
        formula.calculate_internal_energy(test_args.t, test_args.s, test_args.p, 100, test_args.mu, test_args.n)


def test_bad_chemical_potential(test_args: Args) -> None:
    mub = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        formula.calculate_internal_energy(test_args.t, test_args.s, test_args.p, test_args.v, mub, test_args.n)
    with raises(TypeError):
        formula.calculate_internal_energy(test_args.t, test_args.s, test_args.p, test_args.v, 100, test_args.n)


def test_bad_particle_count(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        formula.calculate_internal_energy(test_args.t, test_args.s, test_args.p, test_args.v, test_args.mu, nb)
