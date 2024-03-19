from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import specific_resistance_of_coaxial_waveguide as resistance_law

# Description
## Parameters of the coaxial waveguide: the radius of the inner wire is 1.35 millimeters, the radius of the outer conductor
## is 9.0 millimeters, the relative permeability of the dielectric is 1, the specific conductivity of the conductor is 59.5e6 siemens per meter.
## The frequency of signal is 2 * pi * 100e6 radians per second. The specific resistance will be 1 ohm per meter.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", [
    "relative_permeability", "frequency", "specific_conductivity", "outer_radius",
    "inner_radius"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permeability = 1
    frequency = Quantity(2 * pi * 100e6 * (1 / units.second))
    specific_conductivity = Quantity(59.5e6 * (units.siemens / units.meter))
    outer_radius = Quantity(9 * units.millimeter)
    inner_radius = Quantity(1.35 * units.millimeter)

    return Args(relative_permeability=relative_permeability,
        frequency=frequency,
        specific_conductivity=specific_conductivity,
        outer_radius=outer_radius,
        inner_radius=inner_radius)


def test_basic_specific_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_specific_resistance(test_args.relative_permeability,
        test_args.frequency, test_args.specific_conductivity, test_args.outer_radius,
        test_args.inner_radius)
    assert_equal(result, 0.258 * units.ohm / units.meter)


def test_bad_relative_permeability(test_args: Args) -> None:
    relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_specific_resistance(relative_permeability,
            test_args.frequency, test_args.specific_conductivity, test_args.outer_radius,
            test_args.inner_radius)


def test_bad_frequency(test_args: Args) -> None:
    frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_specific_resistance(test_args.relative_permeability,
            frequency, test_args.specific_conductivity, test_args.outer_radius,
            test_args.inner_radius)
    with raises(TypeError):
        resistance_law.calculate_specific_resistance(test_args.relative_permeability, 100,
            test_args.specific_conductivity, test_args.outer_radius, test_args.inner_radius)


def test_bad_specific_conductivity(test_args: Args) -> None:
    specific_conductivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_specific_resistance(test_args.relative_permeability,
            test_args.frequency, specific_conductivity, test_args.outer_radius,
            test_args.inner_radius)
    with raises(TypeError):
        resistance_law.calculate_specific_resistance(test_args.relative_permeability,
            test_args.frequency, 100, test_args.outer_radius,
            test_args.inner_radius)


def test_bad_radius(test_args: Args) -> None:
    bad_radius = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_specific_resistance(test_args.relative_permeability, test_args.frequency, test_args.specific_conductivity, bad_radius, test_args.inner_radius)
    with raises(TypeError):
        resistance_law.calculate_specific_resistance(test_args.relative_permeability, test_args.frequency, test_args.specific_conductivity, 100, test_args.inner_radius)
    with raises(errors.UnitsError):
        resistance_law.calculate_specific_resistance(test_args.relative_permeability, test_args.frequency, test_args.specific_conductivity, test_args.outer_radius, bad_radius)
    with raises(TypeError):
        resistance_law.calculate_specific_resistance(test_args.relative_permeability, test_args.frequency, test_args.specific_conductivity, test_args.outer_radius, 100)
    with raises(ValueError):
        resistance_law.calculate_specific_resistance(test_args.relative_permeability, test_args.frequency, test_args.specific_conductivity, test_args.inner_radius, test_args.outer_radius)
