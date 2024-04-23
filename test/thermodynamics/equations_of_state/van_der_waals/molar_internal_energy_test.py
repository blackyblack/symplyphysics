from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    prefixes,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.equations_of_state.van_der_waals import molar_internal_energy as energy_law

# Description
## For water in the model of van der Waals equation of state, the molar internal energy amounts
## to -8.43 kJ/mol at temperature T = 300 K. The isochoric molar heat capacity is 74.5 J/(mol*K),
## the bonding forces parameter of van der Waals equation is 5.54 bar*(L/mol)**2 and the molar
## volume of water is 18 cm**3/mol.

Args = namedtuple("Args", "cv t a vm")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cv = Quantity(74.5 * units.joule / (units.mole * units.kelvin))
    t = Quantity(300 * units.kelvin)
    a = Quantity(5.54 * units.bar * (units.liter / units.mole)**2)
    vm = Quantity(18 * units.centimeter**3 / units.mole)
    return Args(cv=cv, t=t, a=a, vm=vm)


def test_law(test_args: Args) -> None:
    result = energy_law.calculate_internal_energy(test_args.cv, test_args.t, test_args.a,
        test_args.vm)
    assert_equal(result, -8.43 * prefixes.kilo * units.joule / units.mole)


def test_bad_molar_heat_capacity(test_args: Args) -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_internal_energy(cb, test_args.t, test_args.a, test_args.vm)
    with raises(TypeError):
        energy_law.calculate_internal_energy(100, test_args.t, test_args.a, test_args.vm)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_internal_energy(test_args.cv, tb, test_args.a, test_args.vm)
    with raises(TypeError):
        energy_law.calculate_internal_energy(test_args.cv, 100, test_args.a, test_args.vm)


def test_bad_bonding_forces_parameter(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_internal_energy(test_args.cv, test_args.t, ab, test_args.vm)
    with raises(TypeError):
        energy_law.calculate_internal_energy(test_args.cv, test_args.t, 100, test_args.vm)


def test_bad_molar_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_internal_energy(test_args.cv, test_args.t, test_args.a, vb)
    with raises(TypeError):
        energy_law.calculate_internal_energy(test_args.cv, test_args.t, test_args.a, 100)
