from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import gibbs_energy_differential as gibbs_law

# Description
## The temperature of a system increases by 1 K and the pressure drops by 10 Pa.
## The entropy of the system is 100 J/K and its volume is 2 m**3. The number of
## particles in the system increases by 2. The chemical potential of the system
## is -10 J. The total change in Gibbs energy amounts to -100 J.

Args = namedtuple("Args", "s dt v dp mu dn")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    s = Quantity(100 * units.joule / units.kelvin)
    dt = Quantity(1 * units.kelvin)
    v = Quantity(2 * units.meter**3)
    dp = Quantity(10 * units.pascal)
    mu = Quantity(-10 * units.joule)
    dn = 2
    return Args(s=s, dt=dt, v=v, dp=dp, mu=mu, dn=dn)


def test_law(test_args: Args) -> None:
    result = gibbs_law.calculate_gibbs_energy_change(test_args.s, test_args.dt, test_args.v, test_args.dp, test_args.mu, test_args.dn)
    assert_equal(result, -100 * units.joule)


def test_bad_entropy(test_args: Args) -> None:
    sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_law.calculate_gibbs_energy_change(sb, test_args.dt, test_args.v, test_args.dp, test_args.mu, test_args.dn)
    with raises(TypeError):
        gibbs_law.calculate_gibbs_energy_change(100, test_args.dt, test_args.v, test_args.dp, test_args.mu, test_args.dn)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_law.calculate_gibbs_energy_change(test_args.s, tb, test_args.v, test_args.dp, test_args.mu, test_args.dn)
    with raises(TypeError):
        gibbs_law.calculate_gibbs_energy_change(test_args.s, 100, test_args.v, test_args.dp, test_args.mu, test_args.dn)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_law.calculate_gibbs_energy_change(test_args.s, test_args.dt, vb, test_args.dp, test_args.mu, test_args.dn)
    with raises(TypeError):
        gibbs_law.calculate_gibbs_energy_change(test_args.s, test_args.dt, 100, test_args.dp, test_args.mu, test_args.dn)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_law.calculate_gibbs_energy_change(test_args.s, test_args.dt, test_args.v, pb, test_args.mu, test_args.dn)
    with raises(TypeError):
        gibbs_law.calculate_gibbs_energy_change(test_args.s, test_args.dt, test_args.v, 100, test_args.mu, test_args.dn)


def test_bad_chemical_potential(test_args: Args) -> None:
    mub = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_law.calculate_gibbs_energy_change(test_args.s, test_args.dt, test_args.v, test_args.dp, mub, test_args.dn)
    with raises(TypeError):
        gibbs_law.calculate_gibbs_energy_change(test_args.s, test_args.dt, test_args.v, test_args.dp, 100, test_args.dn)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_law.calculate_gibbs_energy_change(test_args.s, test_args.dt, test_args.v, test_args.dp, test_args.mu, nb)
