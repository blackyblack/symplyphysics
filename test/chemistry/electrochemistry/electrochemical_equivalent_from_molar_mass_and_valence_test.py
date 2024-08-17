from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.chemistry.electrochemistry import electrochemical_equivalent_from_molar_mass_and_valence as equivalent_law

# Description
## Consider molybdenum with a molar mass of 96 gram per mol. With a valence of 4,
## the electrochemical equivalent is 248 microgram per coulomb.
## https://www.indigomath.ru//raschety/ZBpXGP.html

Args = namedtuple("Args", ["molar_mass", "valence"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    molar_mass = Quantity(96 * (units.gram / units.mol))
    valence = 4
    return Args(molar_mass=molar_mass, valence=valence)


def test_basic_equivalent(test_args: Args) -> None:
    result = equivalent_law.calculate_equivalent(test_args.molar_mass, test_args.valence)
    assert_equal(result, 248.74 * prefixes.micro * units.gram / units.coulomb)


def test_bad_molar_mass(test_args: Args) -> None:
    molar_mass = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        equivalent_law.calculate_equivalent(molar_mass, test_args.valence)
    with raises(TypeError):
        equivalent_law.calculate_equivalent(100, test_args.valence)


def test_bad_valence(test_args: Args) -> None:
    valence = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        equivalent_law.calculate_equivalent(test_args.molar_mass, valence)
    with raises(TypeError):
        equivalent_law.calculate_equivalent(test_args.molar_mass, True)
    with raises(ValueError):
        equivalent_law.calculate_equivalent(test_args.molar_mass, 4.1)
    with raises(ValueError):
        equivalent_law.calculate_equivalent(test_args.molar_mass, -4)
