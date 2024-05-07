from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
    prefixes
)

from symplyphysics.laws.electricity.transmission_lines.microstrip_lines import inductance_of_the_microstrip_line_strip as inductance_law

## The strip width of the microstrip line is 1 millimeter, and the strip length is 5 millimeters, the strip
## thickness is 50 micrometer. Then the inductance of the microstrip line will be equal to 14 microhenry.

Args = namedtuple("Args", ["strip_thickness", "strip_length", "strip_width"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    strip_thickness = Quantity(50 * units.micrometer)
    strip_length = Quantity(5 * units.millimeter)
    strip_width = Quantity(1 * units.millimeter)
    return Args(strip_thickness=strip_thickness,
        strip_length=strip_length,
        strip_width=strip_width)


def test_basic_inductance(test_args: Args) -> None:
    result = inductance_law.calculate_inductance(test_args.strip_thickness,
        test_args.strip_length, test_args.strip_width)
    assert_equal(result, 14 * prefixes.micro * units.henry)


def test_bad_strip_thickness(test_args: Args) -> None:
    bad_strip_thickness = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inductance_law.calculate_inductance(bad_strip_thickness, test_args.strip_length,
            test_args.strip_width)
    with raises(TypeError):
        inductance_law.calculate_inductance(100, test_args.strip_length,
            test_args.strip_width)


def test_bad_strip_length(test_args: Args) -> None:
    bad_strip_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inductance_law.calculate_inductance(test_args.strip_thickness, bad_strip_length,
            test_args.strip_width)
    with raises(TypeError):
        inductance_law.calculate_inductance(test_args.strip_thickness, 100,
            test_args.strip_width)


def test_bad_strip_width(test_args: Args) -> None:
    bad_strip_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inductance_law.calculate_inductance(test_args.strip_thickness, test_args.strip_length,
            bad_strip_width)
    with raises(TypeError):
        inductance_law.calculate_inductance(test_args.strip_thickness, test_args.strip_length,
            100)
