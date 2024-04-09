from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    prefixes,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    internal_energy_change_via_entropy_change_and_volume_change as internal_energy_law,
)

# Description
## An ensamble of particles in a closed reservoir is in thermodynamic equilibrium with the environment.
## The temperature is 300 K and the pressure in the reservoir is 0.1 MPa. The entropy of the ensamble
## has dropped by 1 J/K and the volume has become bigger by 1 mm**3. With chemical potential of the system
## equal to 20 J, the particle count increased by 1. Then the change in internal energy of the system 
## amounts to -280 J.

Args = namedtuple("Args", "t ds p dv mu dn")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(300 * units.kelvin)
    ds = Quantity(-1 * units.joule / units.kelvin)
    p = Quantity(0.1 * prefixes.mega * units.pascal)
    dv = Quantity(1 * units.millimeter**3)
    mu = Quantity(20 * units.joule)
    dn = 1
    return Args(t=t, ds=ds, p=p, dv=dv, mu=mu, dn=dn)


def test_law(test_args: Args) -> None:
    result = internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds,
        test_args.p, test_args.dv, test_args.mu, test_args.dn)
    assert_equal(result, -280 * units.joule)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        internal_energy_law.calculate_internal_energy_change(tb, test_args.ds, test_args.p,
            test_args.dv, test_args.mu, test_args.dn)
    with raises(TypeError):
        internal_energy_law.calculate_internal_energy_change(100, test_args.ds, test_args.p,
            test_args.dv, test_args.mu, test_args.dn)


def test_bad_entropy(test_args: Args) -> None:
    sb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, sb, test_args.p,
            test_args.dv, test_args.mu, test_args.dn)
    with raises(TypeError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, 100, test_args.p,
            test_args.dv, test_args.mu, test_args.dn)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, pb,
            test_args.dv, test_args.mu, test_args.dn)
    with raises(TypeError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, 100,
            test_args.dv, test_args.mu, test_args.dn)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, test_args.p,
            vb, test_args.mu, test_args.dn)
    with raises(TypeError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, test_args.p,
            100, test_args.mu, test_args.dn)


def test_bad_chemical_potential(test_args: Args) -> None:
    mub = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, test_args.p,
            test_args.dv, mub, test_args.dn)
    with raises(TypeError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, test_args.p,
            test_args.dv, 100, test_args.dn)


def test_bad_particle_count(test_args: Args) -> None:
    nb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, test_args.p,
            test_args.dv, test_args.mu, nb)
