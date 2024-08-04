from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.circuits.transmission_lines.microstrip_lines import attenuation_coefficient_in_metal_of_microstrip_line_for_width_greater_thickness as coefficient_law

## The surface resistance of the metal strip is 2.576 milliohm, the effective permittivity is 2.5. The wave resistance of the line is 70 ohm. The thickness
## of the substrate and the metal strip is 7 millimeter and 50 micrometer, respectively. The effective width of the line is 10 millimeter.
## Then the losses in the metal strip will be equal to 0.0226 [1 / meter].

Args = namedtuple("Args", [
    "surface_resistance", "wave_resistance", "thickness_of_substrate", "effective_width",
    "strip_thickness", "effective_permittivity"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    surface_resistance = Quantity(2.576 * prefixes.milli * units.ohm)
    wave_resistance = Quantity(70 * units.ohm)
    thickness_of_substrate = Quantity(7 * units.millimeter)
    effective_width = Quantity(10 * units.millimeter)
    strip_thickness = Quantity(50 * units.micrometer)
    effective_permittivity = 2.5
    return Args(surface_resistance=surface_resistance,
        wave_resistance=wave_resistance,
        thickness_of_substrate=thickness_of_substrate,
        effective_width=effective_width,
        strip_thickness=strip_thickness,
        effective_permittivity=effective_permittivity)


def test_basic_attenuation_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
        test_args.wave_resistance, test_args.thickness_of_substrate, test_args.effective_width,
        test_args.strip_thickness, test_args.effective_permittivity)
    assert_equal(result, 0.0226 * (1 / units.meter))


def test_bad_resistances(test_args: Args) -> None:
    bad_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(bad_resistance, test_args.wave_resistance,
            test_args.thickness_of_substrate, test_args.effective_width, test_args.strip_thickness,
            test_args.effective_permittivity)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(100, test_args.wave_resistance,
            test_args.thickness_of_substrate, test_args.effective_width, test_args.strip_thickness,
            test_args.effective_permittivity)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            bad_resistance, test_args.thickness_of_substrate, test_args.effective_width,
            test_args.strip_thickness, test_args.effective_permittivity)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance, 100,
            test_args.thickness_of_substrate, test_args.effective_width, test_args.strip_thickness,
            test_args.effective_permittivity)


def test_bad_thickness_of_substrate(test_args: Args) -> None:
    thickness_of_substrate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.wave_resistance, thickness_of_substrate, test_args.effective_width,
            test_args.strip_thickness, test_args.effective_permittivity)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.wave_resistance, 100, test_args.effective_width, test_args.strip_thickness,
            test_args.effective_permittivity)
    with raises(ValueError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.wave_resistance, test_args.effective_width, test_args.thickness_of_substrate,
            test_args.strip_thickness, test_args.effective_permittivity)


def test_bad_effective_width(test_args: Args) -> None:
    effective_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.wave_resistance, test_args.thickness_of_substrate, effective_width,
            test_args.strip_thickness, test_args.effective_permittivity)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.wave_resistance, test_args.thickness_of_substrate, 100,
            test_args.strip_thickness, test_args.effective_permittivity)


def test_bad_strip_thickness(test_args: Args) -> None:
    strip_thickness = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.wave_resistance, test_args.thickness_of_substrate, test_args.effective_width,
            strip_thickness, test_args.effective_permittivity)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.wave_resistance, test_args.thickness_of_substrate, test_args.effective_width,
            100, test_args.effective_permittivity)


def test_bad_effective_permittivity(test_args: Args) -> None:
    effective_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.surface_resistance,
            test_args.wave_resistance, test_args.thickness_of_substrate, test_args.effective_width,
            test_args.strip_thickness, effective_permittivity)
