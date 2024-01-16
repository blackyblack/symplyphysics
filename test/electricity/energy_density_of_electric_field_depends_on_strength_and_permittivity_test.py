from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.electricity import energy_density_of_electric_field_depends_on_strength_and_permittivity as energy_density_law

# Description
## It is known that with a permittivity equal to 5 and an electric field intensity equal to 10 Volt / meter,
## energy density of the electric field is 2.2e-9 Joule / meter**3.
## https://www.calculatoratoz.com/ru/energy-density-given-electric-field-calculator/Calc-2224


@fixture(name="test_args")
def test_args_fixture():
    relative_permittivity = 5
    electric_intensity = Quantity(10 * (units.volt / units.meter))

    Args = namedtuple("Args", ["relative_permittivity", "electric_intensity"])
    return Args(relative_permittivity=relative_permittivity, electric_intensity=electric_intensity)


def test_basic_energy_density(test_args):
    result = energy_density_law.calculate_energy_density(test_args.relative_permittivity, test_args.electric_intensity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy / units.volume)
    result = convert_to(result, units.joule / units.meter**3).evalf(5)
    assert result == approx(2.2e-9, rel=0.01)


def test_bad_relative_permittivity(test_args):
    relative_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_density_law.calculate_energy_density(relative_permittivity, test_args.electric_intensity)
    with raises(TypeError):
        energy_density_law.calculate_energy_density(True, test_args.electric_intensity)


def test_bad_electric_intensity(test_args):
    electric_intensity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_density_law.calculate_energy_density(test_args.relative_permittivity, electric_intensity)
    with raises(TypeError):
        energy_density_law.calculate_energy_density(test_args.relative_permittivity, 100)
