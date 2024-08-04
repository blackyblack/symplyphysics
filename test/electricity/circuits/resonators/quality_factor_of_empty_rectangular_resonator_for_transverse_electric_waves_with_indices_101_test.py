from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.circuits.resonators import quality_factor_of_empty_rectangular_resonator_for_transverse_electric_waves_with_indices_101 as factor_law

# Description
## The height, width and length of the resonator are 2 centimeter, respectively.
## The frequency is 2 * pi * 100e6 radian per second.
## The surface resistance and magnetic permeability are equal 2.576 milliohm and 1, respectively.
## Then the quality factor is equal to 1021.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", [
    "angular_frequency", "relative_permeability", "surface_resistance", "resonator_width",
    "resonator_height", "resonator_length"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    angular_frequency = Quantity(2 * pi * 100e6 * (units.radian / units.second))
    relative_permeability = 1
    surface_resistance = Quantity(2.576 * prefixes.milli * units.ohm)
    resonator_width = Quantity(2 * units.centimeter)
    resonator_height = Quantity(2 * units.centimeter)
    resonator_length = Quantity(2 * units.centimeter)

    return Args(angular_frequency=angular_frequency,
        relative_permeability=relative_permeability,
        surface_resistance=surface_resistance,
        resonator_width=resonator_width,
        resonator_height=resonator_height,
        resonator_length=resonator_length)


def test_basic_quality_factor(test_args: Args) -> None:
    result = factor_law.calculate_quality_factor(test_args.angular_frequency,
        test_args.relative_permeability, test_args.surface_resistance,
        (test_args.resonator_width, test_args.resonator_height, test_args.resonator_length))
    assert_equal(result, 1021)


def test_bad_angular_frequency(test_args: Args) -> None:
    angular_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(angular_frequency, test_args.relative_permeability,
            test_args.surface_resistance,
            (test_args.resonator_width, test_args.resonator_height, test_args.resonator_length))
    with raises(TypeError):
        factor_law.calculate_quality_factor(100, test_args.relative_permeability,
            test_args.surface_resistance,
            (test_args.resonator_width, test_args.resonator_height, test_args.resonator_length))


def test_bad_relative_permeability(test_args: Args) -> None:
    relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.angular_frequency, relative_permeability,
            test_args.surface_resistance,
            (test_args.resonator_width, test_args.resonator_height, test_args.resonator_length))


def test_bad_surface_resistance(test_args: Args) -> None:
    bad_surface_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.angular_frequency,
            test_args.relative_permeability, bad_surface_resistance,
            (test_args.resonator_width, test_args.resonator_height, test_args.resonator_length))
    with raises(TypeError):
        factor_law.calculate_quality_factor(test_args.angular_frequency,
            test_args.relative_permeability, 100,
            (test_args.resonator_width, test_args.resonator_height, test_args.resonator_length))


def test_bad_resonator_dimensions(test_args: Args) -> None:
    bad_dimension = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.angular_frequency,
            test_args.relative_permeability, test_args.surface_resistance,
            (bad_dimension, test_args.resonator_height, test_args.resonator_length))
    with raises(TypeError):
        factor_law.calculate_quality_factor(test_args.angular_frequency,
            test_args.relative_permeability, test_args.surface_resistance,
            (100, test_args.resonator_height, test_args.resonator_length))
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.angular_frequency,
            test_args.relative_permeability, test_args.surface_resistance,
            (test_args.resonator_width, bad_dimension, test_args.resonator_length))
    with raises(TypeError):
        factor_law.calculate_quality_factor(test_args.angular_frequency,
            test_args.relative_permeability, test_args.surface_resistance,
            (test_args.resonator_width, 100, test_args.resonator_length))
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.angular_frequency,
            test_args.relative_permeability, test_args.surface_resistance,
            (test_args.resonator_width, test_args.resonator_height, bad_dimension))
    with raises(TypeError):
        factor_law.calculate_quality_factor(test_args.angular_frequency,
            test_args.relative_permeability, test_args.surface_resistance,
            (test_args.resonator_width, test_args.resonator_height, 100))
