from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.chemistry import cross_section_of_interaction_in_recharge_model as cross_section_law

## The ionization energy of atoms is equal to 14.69 electronvolt. The mass of an atom is 4.64e-26 kilogram.
## The pressure is 1 pascal, the temperature is 573 kelvin. The electric field strength is 5857 volt per meter.
## Then the cross-sectional area of the interaction will be 1.31e-18 [meter^2].

Args = namedtuple("Args", [
    "ionization_energy", "mass_of_atom", "pressure", "temperature", "electric_intensity"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    ionization_energy = Quantity(14.69 * units.electronvolt)
    mass_of_atom = Quantity(4.64e-26 * units.kilogram)
    pressure = Quantity(1 * units.pascal)
    temperature = Quantity(573 * units.kelvin)
    electric_intensity = Quantity(5857 * (units.volt / units.meter))

    return Args(ionization_energy=ionization_energy,
        mass_of_atom=mass_of_atom,
        pressure=pressure,
        temperature=temperature,
        electric_intensity=electric_intensity)


def test_basic_cross_sectional_area_of_interaction(test_args: Args) -> None:
    result = cross_section_law.calculate_cross_sectional_area_of_interaction(
        test_args.ionization_energy, test_args.mass_of_atom, test_args.pressure,
        test_args.temperature, test_args.electric_intensity,)
    assert_equal(result, 1.31e-18 * units.meter**2, tolerance=0.01)


def test_bad_ionization_energy(test_args: Args) -> None:
    ionization_energy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(ionization_energy,
            test_args.mass_of_atom, test_args.pressure,
            test_args.temperature, test_args.electric_intensity)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(100,
            test_args.mass_of_atom, test_args.pressure,
            test_args.temperature, test_args.electric_intensity)


def test_bad_mass_of_atom(test_args: Args) -> None:
    mass_of_atom = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.ionization_energy,
            mass_of_atom, test_args.pressure, test_args.temperature, test_args.electric_intensity)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.ionization_energy,
            100, test_args.pressure, test_args.temperature, test_args.electric_intensity)


def test_bad_pressure(test_args: Args) -> None:
    pressure = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.ionization_energy,
            test_args.mass_of_atom, pressure, test_args.temperature, test_args.electric_intensity)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.ionization_energy,
            test_args.mass_of_atom, 100, test_args.temperature, test_args.electric_intensity)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.ionization_energy,
            test_args.mass_of_atom, test_args.pressure, temperature, test_args.electric_intensity)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.ionization_energy,
            test_args.mass_of_atom, test_args.pressure, 100, test_args.electric_intensity)


def test_bad_electric_intensity(test_args: Args) -> None:
    electric_intensity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.ionization_energy,
            test_args.mass_of_atom, test_args.pressure, test_args.temperature, electric_intensity)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(test_args.ionization_energy,
            test_args.mass_of_atom, test_args.pressure, test_args.temperature, 100)
