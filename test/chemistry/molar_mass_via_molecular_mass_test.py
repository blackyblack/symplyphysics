from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.chemistry import (
    molar_mass_via_molecular_mass as molar_mass_law,)

# Description
## The molar mass of Argon (m = 39.948 amu) is 39.948 g/mol.

Args = namedtuple("Args", "m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(39.948 * units.amu)
    return Args(m=m)


def test_law(test_args: Args) -> None:
    result = molar_mass_law.calculate_molar_mass(test_args.m)
    assert_equal(result, 39.948 * units.gram / units.mole)


def test_bad_mass() -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        molar_mass_law.calculate_molar_mass(mb)
    with raises(TypeError):
        molar_mass_law.calculate_molar_mass(100)
