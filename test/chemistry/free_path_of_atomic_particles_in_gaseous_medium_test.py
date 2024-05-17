from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import free_path_of_atomic_particles_in_gaseous_medium as path_law

# Description
## The pressure is 1.43 pascal. The temperature is 573 kelvin. The cross-sectional area of the interaction is 2.41e-19 meter^2.
## Then the free path length is 22.96 millimeter.

Args = namedtuple("Args", ["pressure", "temperature", "cross_sectional_area_of_interaction"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    pressure = Quantity(1.43 * units.pascal)
    temperature = Quantity(573 * units.kelvin)
    cross_sectional_area_of_interaction = Quantity(2.41e-19 * (units.meter**2))

    return Args(pressure=pressure,
        temperature=temperature,
        cross_sectional_area_of_interaction=cross_sectional_area_of_interaction)


def test_basic_free_path_length(test_args: Args) -> None:
    result = path_law.calculate_free_path_length(test_args.pressure, test_args.temperature,
        test_args.cross_sectional_area_of_interaction)
    assert_equal(result, 22.96 * units.millimeter)


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
