from collections import namedtuple
from pytest import approx, fixture, raises

from sympy import cos, pi
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.condensed_matter import effective_mass_of_the_electron_from_the_energy as ef_mass_el

# Description
## Let's consider the one-dimensional model for an electron presented on the website:
## https://www.slideserve.com/kiele/lecture-ix-dr-hab-ewa-popko-electron-dynamics
## For values of k such that k << pi/a - k is wavenumber, a is period of the structure - the effective mass
## of the electron is found. The effective mass is 0.24*m_e, where m_e is the usual mass of an electron equal
## to 9.109e-31 kilograms.


@fixture(name="test_args")
def test_args_fixture():
    period_structure_not_quantity = 2e-10
    period_structure = Quantity(period_structure_not_quantity * (units.meter))
    gamma = Quantity(4 * 1.6e-19 * (units.joule))
    wavenumber = Quantity((pi / period_structure_not_quantity / 1e3) * (1/units.meter))
    energy_function = -2 * gamma * cos(wavenumber * period_structure)

    increased_wavenumber = Quantity(wavenumber * 3)
    energy_function_new = -2 * gamma * cos(increased_wavenumber * period_structure)

    mass_electron = Quantity(9.109e-31 * units.kilogram)

    Args = namedtuple("Args", ["energy_function", "energy_function_new",
                               "wavenumber", "increased_wavenumber", "mass_electron"])
    return Args(energy_function=energy_function,
                energy_function_new=energy_function_new,
                wavenumber=wavenumber,
                increased_wavenumber=increased_wavenumber,
                mass_electron=mass_electron)


def test_basic_mass(test_args):
    result = ef_mass_el.calculate_mass(test_args.energy_function, test_args.wavenumber)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass)
    result = convert_to(result, units.kilogram).evalf(6)
    assert result == approx(0.24 * convert_to(test_args.mass_electron, units.kilogram).evalf(6), 0.00001)


def test_mass_increases_with_increasing_wave_number(test_args):
    result_wavenumber = ef_mass_el.calculate_mass(test_args.energy_function, test_args.wavenumber)
    result_increased_wavenumber = ef_mass_el.calculate_mass(test_args.energy_function_new, test_args.increased_wavenumber)
    result_wavenumber = convert_to(result_wavenumber, units.kilogram).evalf(6)
    result_increased_wavenumber = convert_to(result_increased_wavenumber, units.kilogram).evalf(6)
    assert result_increased_wavenumber - result_wavenumber > 0


def test_bad_propagation_vec(test_args):
    wavenumber = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ef_mass_el.calculate_mass(test_args.energy_function, wavenumber)
    with raises(TypeError):
        ef_mass_el.calculate_mass(test_args.energy_function, 100)
