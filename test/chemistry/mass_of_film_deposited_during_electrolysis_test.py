from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.chemistry import mass_of_film_deposited_during_electrolysis as mass_law

# Description
## The current is 10 milliampere. The molar mass of a metal atom is 63.546 gram per mole. The current output is 1.05, and the valence of the metal atom is 29.
## The time is 1800 seconds. Then the mass of the resulting film will be equal to 0.429 milligram.

Args = namedtuple("Args", ["current", "molar_mass", "current_output", "valence", "time"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    current = Quantity(10 * prefixes.milli * units.ampere)
    molar_mass = Quantity(63.546 * units.gram / units.mol)
    current_output = 1.05
    valence = 29
    time = Quantity(1800 * units.second)

    return Args(current=current,
                molar_mass=molar_mass,
                current_output=current_output,
                valence=valence,
                time=time)


def test_basic_mass_of_film(test_args: Args) -> None:
    result = mass_law.calculate_mass_of_film(test_args.current, test_args.molar_mass, test_args.current_output, test_args.valence, test_args.time)
    assert_equal(result, 0.429 * units.milligram)


def test_bad_current(test_args: Args) -> None:
    current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mass_law.calculate_mass_of_film(current, test_args.molar_mass, test_args.current_output, test_args.valence, test_args.time)
    with raises(TypeError):
        mass_law.calculate_mass_of_film(100, test_args.molar_mass, test_args.current_output, test_args.valence, test_args.time)


def test_bad_molar_mass(test_args: Args) -> None:
    molar_mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mass_law.calculate_mass_of_film(test_args.current, molar_mass, test_args.current_output, test_args.valence, test_args.time)
    with raises(TypeError):
        mass_law.calculate_mass_of_film(test_args.current, 100, test_args.current_output, test_args.valence, test_args.time)


def test_bad_current_output(test_args: Args) -> None:
    current_output = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mass_law.calculate_mass_of_film(test_args.current, test_args.molar_mass, current_output, test_args.valence, test_args.time)


def test_bad_valence(test_args: Args) -> None:
    valence = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mass_law.calculate_mass_of_film(test_args.current, test_args.molar_mass, test_args.current_output, valence, test_args.time)


def test_bad_time(test_args: Args) -> None:
    time = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mass_law.calculate_mass_of_film(test_args.current, test_args.molar_mass, test_args.current_output, test_args.valence, time)
    with raises(TypeError):
        mass_law.calculate_mass_of_film(test_args.current, test_args.molar_mass, test_args.current_output, test_args.valence, 100)
