from collections import namedtuple
from pytest import fixture, raises
from sympy import cos, pi
from symplyphysics import (
    assert_equal,
    errors,
    units,
    convert_to,
    Quantity,
)
from symplyphysics.laws.condensed_matter import effective_mass_of_electron_via_energy as ef_mass_el

# Description
## Let's consider the one-dimensional model for an electron presented on the website:
## https://www.slideserve.com/kiele/lecture-ix-dr-hab-ewa-popko-electron-dynamics
## For values of k such that k << pi/a - k is wavenumber, a is period of the structure - the effective mass
## of the electron is found. The effective mass is 0.24*m_e, where m_e is the usual mass of an electron equal
## to 9.109e-31 kilograms.

Args = namedtuple("Args", [
    "energy_function", "energy_function_new", "wavenumber", "increased_wavenumber", "mass_electron"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    period_structure = Quantity(2e-10 * units.meter)
    gamma = Quantity(4 * 1.6e-19 * units.joule)
    wavenumber = Quantity(pi / period_structure / 1000)
    energy_function = -2 * gamma * cos(wavenumber * period_structure)
    increased_wavenumber = Quantity(wavenumber * 3)
    energy_function_new = -2 * gamma * cos(increased_wavenumber * period_structure)
    mass_electron = Quantity(9.109e-31 * units.kilogram)

    return Args(energy_function=energy_function,
        energy_function_new=energy_function_new,
        wavenumber=wavenumber,
        increased_wavenumber=increased_wavenumber,
        mass_electron=mass_electron)


def test_basic_mass(test_args: Args) -> None:
    result = ef_mass_el.calculate_mass(test_args.energy_function, test_args.wavenumber)
    assert_equal(result, 0.24 * test_args.mass_electron, tolerance=0.01)


def test_mass_increases_with_increasing_wave_number(test_args: Args) -> None:
    result_wavenumber = ef_mass_el.calculate_mass(test_args.energy_function, test_args.wavenumber)
    result_increased_wavenumber = ef_mass_el.calculate_mass(test_args.energy_function_new,
        test_args.increased_wavenumber)
    result_wavenumber = convert_to(result_wavenumber, units.kilogram).evalf(6)
    result_increased_wavenumber = convert_to(result_increased_wavenumber, units.kilogram).evalf(6)
    assert result_increased_wavenumber - result_wavenumber > 0


def test_bad_propagation_vec(test_args: Args) -> None:
    wavenumber = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ef_mass_el.calculate_mass(test_args.energy_function, wavenumber)
    with raises(TypeError):
        ef_mass_el.calculate_mass(test_args.energy_function, 100)
