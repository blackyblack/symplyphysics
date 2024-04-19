from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    free_energy_differential as free_energy_law,)

# Description
## A closed system in thermal equlibrium has entropy S = 1 J/k and pressure p = 1 Pa inside.
## The temperature of the system drops by 0.01 K and the volume grows by 0.02 m**3. The number
## of particles decreases by 1, the chemical potential is 0.1 J.
## Then the change in Helmholtz free energy of the system is -0.11 J.

Args = namedtuple("Args", "s dt p dv mu dn")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    s = Quantity(1 * units.joule / units.kelvin)
    dt = Quantity(-0.01 * units.kelvin)
    p = Quantity(1 * units.pascal)
    dv = Quantity(0.02 * units.meter**3)
    mu = Quantity(0.1 * units.joule)
    dn = -1
    return Args(s=s, dt=dt, p=p, dv=dv, mu=mu, dn=dn)


def test_law(test_args: Args) -> None:
    result = free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, test_args.p,
        test_args.dv, test_args.mu, test_args.dn)
    assert_equal(result, -0.11 * units.joule)


def test_bad_entropy(test_args: Args) -> None:
    sb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        free_energy_law.calculate_free_energy_change(sb, test_args.dt, test_args.p, test_args.dv,
            test_args.mu, test_args.dn)
    with raises(TypeError):
        free_energy_law.calculate_free_energy_change(100, test_args.dt, test_args.p, test_args.dv,
            test_args.mu, test_args.dn)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        free_energy_law.calculate_free_energy_change(test_args.s, tb, test_args.p, test_args.dv,
            test_args.mu, test_args.dn)
    with raises(TypeError):
        free_energy_law.calculate_free_energy_change(test_args.s, 100, test_args.p, test_args.dv,
            test_args.mu, test_args.dn)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, pb, test_args.dv,
            test_args.mu, test_args.dn)
    with raises(TypeError):
        free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, 100, test_args.dv,
            test_args.mu, test_args.dn)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, test_args.p, vb,
            test_args.mu, test_args.dn)
    with raises(TypeError):
        free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, test_args.p, 100,
            test_args.mu, test_args.dn)


def test_bad_chemical_potential(test_args: Args) -> None:
    mub = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, test_args.p,
            test_args.dv, mub, test_args.dn)
    with raises(TypeError):
        free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, test_args.p,
            test_args.dv, 100, test_args.dn)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, test_args.p,
            test_args.dv, test_args.mu, nb)
