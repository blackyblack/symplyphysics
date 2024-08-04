from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits.waveguides import power_carried_by_coaxial_waveguide as power_law

# Description
## Parameters of the coaxial waveguide: the diameter of the inner wire is 1.35 millimeters, the diameter of the outer conductor
## is 9.0 millimeters, the relative permittivity of the dielectric is 2.2, the relative permeability of the dielectric is 1.
## The voltage is 10 volts. Then the transferred power will be equal to 0.897 watt.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", [
    "relative_permittivity", "relative_permeability", "voltage", "outer_diameter", "inner_diameter"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 2.2
    relative_permeability = 1
    voltage = Quantity(10 * units.volt)
    outer_diameter = Quantity(9 * units.millimeter)
    inner_diameter = Quantity(1.35 * units.millimeter)

    return Args(relative_permittivity=relative_permittivity,
        relative_permeability=relative_permeability,
        voltage=voltage,
        outer_diameter=outer_diameter,
        inner_diameter=inner_diameter)


def test_basic_waveguide_power(test_args: Args) -> None:
    result = power_law.calculate_waveguide_power(test_args.relative_permittivity,
        test_args.relative_permeability, test_args.voltage, test_args.outer_diameter,
        test_args.inner_diameter)
    assert_equal(result, 0.897 * units.watt)


def test_bad_relative_permittivity(test_args: Args) -> None:
    relative_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_waveguide_power(relative_permittivity, test_args.relative_permeability,
            test_args.voltage, test_args.outer_diameter, test_args.inner_diameter)


def test_bad_relative_permeability(test_args: Args) -> None:
    relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_waveguide_power(test_args.relative_permittivity, relative_permeability,
            test_args.voltage, test_args.outer_diameter, test_args.inner_diameter)


def test_bad_voltage(test_args: Args) -> None:
    voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_waveguide_power(test_args.relative_permittivity,
            test_args.relative_permeability, voltage, test_args.outer_diameter,
            test_args.inner_diameter)
    with raises(TypeError):
        power_law.calculate_waveguide_power(test_args.relative_permittivity,
            test_args.relative_permeability, 100, test_args.outer_diameter,
            test_args.inner_diameter)


def test_bad_diameter(test_args: Args) -> None:
    bad_diameter = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_waveguide_power(test_args.relative_permittivity,
            test_args.relative_permeability, test_args.voltage, bad_diameter,
            test_args.inner_diameter)
    with raises(TypeError):
        power_law.calculate_waveguide_power(test_args.relative_permittivity,
            test_args.relative_permeability, test_args.voltage, 100, test_args.inner_diameter)
    with raises(errors.UnitsError):
        power_law.calculate_waveguide_power(test_args.relative_permittivity,
            test_args.relative_permeability, test_args.voltage, test_args.outer_diameter,
            bad_diameter)
    with raises(TypeError):
        power_law.calculate_waveguide_power(test_args.relative_permittivity,
            test_args.relative_permeability, test_args.voltage, test_args.outer_diameter, 100)
    with raises(ValueError):
        power_law.calculate_waveguide_power(test_args.relative_permittivity,
            test_args.relative_permeability, test_args.voltage, test_args.inner_diameter,
            test_args.outer_diameter)
