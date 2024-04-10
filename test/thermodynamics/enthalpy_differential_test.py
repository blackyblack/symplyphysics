from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.thermodynamics import (
    enthalpy_differential as enthalpy_law,)

# Description
## A closed homogeneous system is in thermal equilibrium. Its temperature is 400 K and
## volume is 0.5 m**3. The system undergoes a reversible process in which its entropy
## drops by 1e-5 J/K and its pressure rises by 0.01 Pa. The number of particles increases
## by 2 and the chemical potential of the system is 1 mJ. The change in system's enthalpy
## amounts to 3 mJ.

Args = namedtuple("Args", "t ds v dp mu dn")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(400 * units.kelvin)
    ds = Quantity(-1e-5 * units.joule / units.kelvin)
    v = Quantity(0.5 * units.meter**3)
    dp = Quantity(0.01 * units.pascal)
    mu = Quantity(1 * prefixes.milli * units.joule)
    dn = 2
    return Args(t=t, ds=ds, v=v, dp=dp, mu=mu, dn=dn)


def test_law(test_args: Args) -> None:
    result = enthalpy_law.calculate_enthalpy_change(test_args.t, test_args.ds, test_args.v,
        test_args.dp, test_args.mu, test_args.dn)
    assert_equal(result, 3 * prefixes.milli * units.joule)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy_change(tb, test_args.ds, test_args.v, test_args.dp, test_args.mu, test_args.dn)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy_change(100, test_args.ds, test_args.v, test_args.dp, test_args.mu, test_args.dn)


def test_bad_entropy(test_args: Args) -> None:
    sb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy_change(test_args.t, sb, test_args.v, test_args.dp, test_args.mu, test_args.dn)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy_change(test_args.t, 100, test_args.v, test_args.dp, test_args.mu, test_args.dn)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy_change(test_args.t, test_args.ds, vb, test_args.dp, test_args.mu, test_args.dn)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy_change(test_args.t, test_args.ds, 100, test_args.dp, test_args.mu, test_args.dn)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy_change(test_args.t, test_args.ds, test_args.v, pb, test_args.mu, test_args.dn)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy_change(test_args.t, test_args.ds, test_args.v, 100, test_args.mu, test_args.dn)


def test_bad_chemical_potential(test_args: Args) -> None:
    mub = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy_change(test_args.t, test_args.ds, test_args.v, test_args.dp, mub, test_args.dn)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy_change(test_args.t, test_args.ds, test_args.v, test_args.dp, 100, test_args.dn)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy_change(test_args.t, test_args.ds, test_args.v, test_args.dp, test_args.mu, nb)
