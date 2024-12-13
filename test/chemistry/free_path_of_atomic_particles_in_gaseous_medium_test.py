from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, quantities)
from symplyphysics.laws.chemistry import free_path_of_atomic_particles_in_gaseous_medium as path_law

# Description
## Under normal conditions (T = 0°C, p = 1 atm) the mean free path of hydrogen molecules is 130 nm.
## The kinetic diameter of H2 molecules is 289 pm, therefore the cross sectional area of interactions
## amounts to 2.62e-19 m**2.

# Links
## Mean free path - https://physicstasks.eu/3943/mean-free-path-of-hydrogen
## Kinetic diameter - https://en.wikipedia.org/wiki/Kinetic_diameter

Args = namedtuple("Args", ["pressure", "temperature", "cross_sectional_area_of_interaction"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    pressure = Quantity(1 * units.atmosphere)
    temperature = quantities.standard_conditions_temperature
    cross_sectional_area_of_interaction = Quantity(2.62e-19 * (units.meter**2))

    return Args(pressure=pressure,
        temperature=temperature,
        cross_sectional_area_of_interaction=cross_sectional_area_of_interaction)


def test_basic_free_path_length(test_args: Args) -> None:
    result = path_law.calculate_free_path_length(test_args.pressure, test_args.temperature,
        test_args.cross_sectional_area_of_interaction)
    assert_equal(result, 130 * units.nanometer, tolerance=0.3)


def test_bad_pressure(test_args: Args) -> None:
    pressure = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        path_law.calculate_free_path_length(pressure, test_args.temperature,
            test_args.cross_sectional_area_of_interaction)
    with raises(TypeError):
        path_law.calculate_free_path_length(100, test_args.temperature,
            test_args.cross_sectional_area_of_interaction)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        path_law.calculate_free_path_length(test_args.pressure, temperature,
            test_args.cross_sectional_area_of_interaction)
    with raises(TypeError):
        path_law.calculate_free_path_length(test_args.pressure, 100,
            test_args.cross_sectional_area_of_interaction)


def test_bad_cross_sectional_area_of_interaction(test_args: Args) -> None:
    cross_sectional_area_of_interaction = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        path_law.calculate_free_path_length(test_args.pressure, test_args.temperature,
            cross_sectional_area_of_interaction)
    with raises(TypeError):
        path_law.calculate_free_path_length(test_args.pressure, test_args.temperature, 100)
