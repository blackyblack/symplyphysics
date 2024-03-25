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
## has dropped by 1 J/K and the volume has become bigger by 1 mm**3. Then the change in internal energy
## of the system amounts to -300 J.

Args = namedtuple("Args", "t ds p dv")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(300 * units.kelvin)
    ds = Quantity(-1 * units.joule / units.kelvin)
    p = Quantity(0.1 * prefixes.mega * units.pascal)
    dv = Quantity(1 * units.millimeter**3)
    return Args(t=t, ds=ds, p=p, dv=dv)


def test_law(test_args: Args) -> Args:
    result = internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, test_args.p, test_args.dv)
    assert_equal(result, -300 * units.joule)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        internal_energy_law.calculate_internal_energy_change(tb, test_args.ds, test_args.p, test_args.dv)
    with raises(TypeError):
        internal_energy_law.calculate_internal_energy_change(100, test_args.ds, test_args.p, test_args.dv)


def test_bad_entropy(test_args: Args) -> None:
    sb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, sb, test_args.p, test_args.dv)
    with raises(TypeError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, 100, test_args.p, test_args.dv)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, pb, test_args.dv)
    with raises(TypeError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, 100, test_args.dv)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, test_args.p, vb)
    with raises(TypeError):
        internal_energy_law.calculate_internal_energy_change(test_args.t, test_args.ds, test_args.p, 100)
