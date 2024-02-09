from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal
)

from symplyphysics.laws.thermodynamics import change_in_entropy_of_ideal_gas_from_volume_and_temperature as change_in_entropy

## Source of numbers: https://edu.tltsu.ru/er/book_view.php?book_id=48e&page_id=3918

@fixture(name="test_args")
def test_args_fixture():
    mass = Quantity(0.008 * units.kilograms)
    molar_mass = Quantity(0.032 * units.kilograms / units.mole)
    molar_heat_capacity = Quantity(29.4 * units.joule / (units.kelvin * units.mole))
    final_temperature = Quantity(573 * units.kelvin)
    initial_temperature = Quantity(353 * units.kelvin)
    final_volume = Quantity(0.04 * units.meter**3)
    start_volume = Quantity(0.01 * units.meter**3)
    Args = namedtuple("Args", ["mass", "molar_mass", "molar_heat_capacity", "final_temperature", "initial_temperature", "final_volume", "start_volume"])
    return Args(
        mass=mass,
        molar_mass=molar_mass,
        molar_heat_capacity=molar_heat_capacity,
        final_temperature=final_temperature,
        initial_temperature=initial_temperature,
        final_volume=final_volume,
        start_volume=start_volume
    )


def test_basic_law(test_args):
    result = change_in_entropy.calculate_entropy_change(test_args.mass, test_args.molar_mass, test_args.molar_heat_capacity, test_args.final_temperature, test_args.initial_temperature, test_args.final_volume, test_args.start_volume)
    assert_equal(result, 6.442 * (units.joule / units.kelvin))


def test_bad_mass(test_args):
    bad_mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        change_in_entropy.calculate_entropy_change(bad_mass, test_args.molar_mass, test_args.molar_heat_capacity, test_args.final_temperature, test_args.initial_temperature, test_args.final_volume, test_args.start_volume)
    with raises(TypeError):
        change_in_entropy.calculate_entropy_change(100, test_args.molar_mass, test_args.molar_heat_capacity, test_args.final_temperature, test_args.initial_temperature, test_args.final_volume, test_args.start_volume)


def test_bad_molar_mass(test_args):
    bad_molar_mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        change_in_entropy.calculate_entropy_change(test_args.mass, bad_molar_mass, test_args.molar_heat_capacity, test_args.final_temperature, test_args.initial_temperature, test_args.final_volume, test_args.start_volume)
    with raises(TypeError):
        change_in_entropy.calculate_entropy_change(test_args.mass, 100, test_args.molar_heat_capacity, test_args.final_temperature, test_args.initial_temperature, test_args.final_volume, test_args.start_volume)


def test_bad_molar_heat_capacity(test_args):
    bad_molar_heat_capacity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        change_in_entropy.calculate_entropy_change(test_args.mass, test_args.molar_mass, bad_molar_heat_capacity, test_args.final_temperature, test_args.initial_temperature, test_args.final_volume, test_args.start_volume)
    with raises(TypeError):
        change_in_entropy.calculate_entropy_change(test_args.mass, test_args.molar_mass, 100, test_args.final_temperature, test_args.initial_temperature, test_args.final_volume, test_args.start_volume)


def test_bad_volume(test_args):
    bad_volume = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        change_in_entropy.calculate_entropy_change(test_args.mass, test_args.molar_mass, test_args.molar_heat_capacity, test_args.final_temperature, test_args.initial_temperature, bad_volume, test_args.start_volume)
    with raises(TypeError):
        change_in_entropy.calculate_entropy_change(test_args.mass, test_args.molar_mass, test_args.molar_heat_capacity, test_args.final_temperature, test_args.initial_temperature, 100, test_args.start_volume)
    with raises(errors.UnitsError):
        change_in_entropy.calculate_entropy_change(test_args.mass, test_args.molar_mass, test_args.molar_heat_capacity, test_args.final_temperature, test_args.initial_temperature, test_args.final_volume, bad_volume)
    with raises(TypeError):
        change_in_entropy.calculate_entropy_change(test_args.mass, test_args.molar_mass, test_args.molar_heat_capacity, test_args.final_temperature, test_args.initial_temperature, test_args.final_volume, 100)


def test_bad_temperature(test_args):
    bad_temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        change_in_entropy.calculate_entropy_change(test_args.mass, test_args.molar_mass, test_args.molar_heat_capacity, bad_temperature, test_args.initial_temperature, test_args.final_volume, test_args.start_volume)
    with raises(TypeError):
        change_in_entropy.calculate_entropy_change(test_args.mass, test_args.molar_mass, test_args.molar_heat_capacity, 100, test_args.initial_temperature, test_args.final_volume, test_args.start_volume)
    with raises(errors.UnitsError):
        change_in_entropy.calculate_entropy_change(test_args.mass, test_args.molar_mass, test_args.molar_heat_capacity, test_args.final_temperature, bad_temperature, test_args.final_volume, test_args.start_volume)
    with raises(TypeError):
        change_in_entropy.calculate_entropy_change(test_args.mass, test_args.molar_mass, test_args.molar_heat_capacity, test_args.final_temperature, 100, test_args.final_volume, test_args.start_volume)
