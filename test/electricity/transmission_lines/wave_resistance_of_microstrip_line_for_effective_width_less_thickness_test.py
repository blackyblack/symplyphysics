from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.transmission_lines import wave_resistance_of_microstrip_line_for_effective_width_less_thickness as resistance_law

## The effective permittivity of the microstrip line is 2.5. The effective width of the microstrip line is 3 millimeter, and the thickness of the substrate
## is 4 millimeters. Then the wave resistance of the microstrip line will be equal to 90.49 ohm.

Args = namedtuple("Args", ["effective_permittivity", "thickness_of_substrate", "effective_width"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    effective_permittivity = 2.5
    thickness_of_substrate = Quantity(4 * units.millimeter)
    effective_width = Quantity(3 * units.millimeter)
    return Args(effective_permittivity=effective_permittivity,
        thickness_of_substrate=thickness_of_substrate,
        effective_width=effective_width)


def test_basic_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_resistance(test_args.effective_permittivity,
        test_args.thickness_of_substrate, test_args.effective_width)
    assert_equal(result, 90.49 * units.ohm)


def test_bad_effective_permittivity(test_args: Args) -> None:
    bad_effective_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(bad_effective_permittivity,
            test_args.thickness_of_substrate, test_args.effective_width)


def test_bad_thickness_of_substrate(test_args: Args) -> None:
    bad_thickness_of_substrate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(test_args.effective_permittivity,
            bad_thickness_of_substrate, test_args.effective_width)
    with raises(TypeError):
        resistance_law.calculate_resistance(test_args.effective_permittivity, 100,
            test_args.effective_width)


def test_bad_effective_width(test_args: Args) -> None:
    bad_effective_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(test_args.effective_permittivity,
            test_args.thickness_of_substrate, bad_effective_width)
    with raises(TypeError):
        resistance_law.calculate_resistance(test_args.effective_permittivity,
            test_args.thickness_of_substrate, 100)
    with raises(ValueError):
        resistance_law.calculate_resistance(test_args.effective_permittivity,
            test_args.effective_width, test_args.thickness_of_substrate)
