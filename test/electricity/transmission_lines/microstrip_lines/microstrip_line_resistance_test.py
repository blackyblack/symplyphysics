from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.transmission_lines.microstrip_lines import microstrip_line_resistance as resistance_law

## The strip width of the microstrip line is 1 millimeter, and the strip length is 5 millimeters, the strip
## thickness is 50 micrometer. The surface resistance of the metal strip is 2.576 milliohm.
## Then the resistance of the microstrip line will be equal to 10.43 milliohm.

Args = namedtuple("Args", ["strip_thickness", "strip_length", "strip_width", "surface_resistance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    strip_thickness = Quantity(50 * units.micrometer)
    strip_length = Quantity(5 * units.millimeter)
    strip_width = Quantity(1 * units.millimeter)
    surface_resistance = Quantity(2.576 * prefixes.milli * units.ohm)
    return Args(strip_thickness=strip_thickness,
        strip_length=strip_length,
        strip_width=strip_width,
        surface_resistance=surface_resistance)


def test_basic_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_resistance(test_args.strip_thickness, test_args.strip_length,
        test_args.strip_width, test_args.surface_resistance)
    assert_equal(result, 10.43 * prefixes.milli * units.ohm)


def test_bad_strip_thickness(test_args: Args) -> None:
    bad_strip_thickness = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(bad_strip_thickness, test_args.strip_length,
            test_args.strip_width, test_args.surface_resistance)
    with raises(TypeError):
        resistance_law.calculate_resistance(100, test_args.strip_length, test_args.strip_width,
            test_args.surface_resistance)


def test_bad_strip_length(test_args: Args) -> None:
    bad_strip_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(test_args.strip_thickness, bad_strip_length,
            test_args.strip_width, test_args.surface_resistance)
    with raises(TypeError):
        resistance_law.calculate_resistance(test_args.strip_thickness, 100, test_args.strip_width,
            test_args.surface_resistance)


def test_bad_strip_width(test_args: Args) -> None:
    bad_strip_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(test_args.strip_thickness, test_args.strip_length,
            bad_strip_width, test_args.surface_resistance)
    with raises(TypeError):
        resistance_law.calculate_resistance(test_args.strip_thickness, test_args.strip_length, 100,
            test_args.surface_resistance)


def test_bad_surface_resistance(test_args: Args) -> None:
    bad_surface_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(test_args.strip_thickness, test_args.strip_length,
            test_args.strip_width, bad_surface_resistance)
    with raises(TypeError):
        resistance_law.calculate_resistance(test_args.strip_thickness, test_args.strip_length,
            test_args.strip_width, 100)
