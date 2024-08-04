from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.circuits.transmission_lines.microstrip_lines import effective_width_for_microstrip_line_for_width_greater_thickness as effective_width_law

## The width of the microstrip line is 1 millimeter, and the thickness of the substrate is 5 millimeters, the strip
## thickness is 50 micrometer. Then the effective width of the microstrip line will be equal to 1.125 millimeter.

Args = namedtuple("Args", ["strip_thickness", "thickness_of_substrate", "width"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    strip_thickness = Quantity(50 * units.micrometer)
    thickness_of_substrate = Quantity(5 * units.millimeter)
    width = Quantity(1 * units.millimeter)
    return Args(strip_thickness=strip_thickness,
        thickness_of_substrate=thickness_of_substrate,
        width=width)


def test_basic_effective_width(test_args: Args) -> None:
    result = effective_width_law.calculate_effective_width(test_args.strip_thickness,
        test_args.thickness_of_substrate, test_args.width)
    assert_equal(result, 1.125 * units.millimeter)


def test_bad_strip_thickness(test_args: Args) -> None:
    bad_strip_thickness = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        effective_width_law.calculate_effective_width(bad_strip_thickness,
            test_args.thickness_of_substrate, test_args.width)


def test_bad_thickness_of_substrate(test_args: Args) -> None:
    bad_thickness_of_substrate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        effective_width_law.calculate_effective_width(test_args.strip_thickness,
            bad_thickness_of_substrate, test_args.width)
    with raises(TypeError):
        effective_width_law.calculate_effective_width(test_args.strip_thickness, 100,
            test_args.width)


def test_bad_width(test_args: Args) -> None:
    bad_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        effective_width_law.calculate_effective_width(test_args.strip_thickness,
            test_args.thickness_of_substrate, bad_width)
    with raises(TypeError):
        effective_width_law.calculate_effective_width(test_args.strip_thickness,
            test_args.thickness_of_substrate, 100)
    thickness_of_substrate = Quantity(7 * units.millimeter)
    with raises(ValueError):
        effective_width_law.calculate_effective_width(test_args.strip_thickness,
            thickness_of_substrate, test_args.width)
