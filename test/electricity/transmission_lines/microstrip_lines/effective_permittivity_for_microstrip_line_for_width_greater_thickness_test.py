from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.transmission_lines.microstrip_lines import effective_permittivity_for_microstrip_line_for_width_greater_thickness as permittivity_law

## The relative permittivity of the microstrip line dielectric is 4. The width of the microstrip line is 5 millimeter, and
## the thickness of the substrate is 4 millimeters, the strip thickness is 50 micrometer.
## Then the effective permittivity is 2.953.

Args = namedtuple("Args", [
    "relative_permittivity", "strip_thickness", "thickness_of_substrate", "width"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 4
    strip_thickness = Quantity(50 * units.micrometer)
    thickness_of_substrate = Quantity(4 * units.millimeter)
    width = Quantity(5 * units.millimeter)
    return Args(relative_permittivity=relative_permittivity,
        strip_thickness=strip_thickness,
        thickness_of_substrate=thickness_of_substrate,
        width=width)


def test_basic_effective_permittivity(test_args: Args) -> None:
    result = permittivity_law.calculate_effective_permittivity(
        test_args.relative_permittivity, test_args.strip_thickness, test_args.thickness_of_substrate,
        test_args.width)
    assert_equal(result, 2.953)


def test_bad_permittivity(test_args: Args) -> None:
    bad_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(bad_permittivity,
            test_args.strip_thickness, test_args.thickness_of_substrate,
            test_args.width)


def test_bad_strip_thickness(test_args: Args) -> None:
    bad_strip_thickness = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            bad_strip_thickness, test_args.thickness_of_substrate, test_args.width)
    with raises(TypeError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            100, test_args.thickness_of_substrate, test_args.width)


def test_bad_thickness_of_substrate(test_args: Args) -> None:
    bad_thickness_of_substrate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            test_args.strip_thickness, bad_thickness_of_substrate, test_args.width)
    with raises(TypeError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            test_args.strip_thickness, 100, test_args.width)


def test_bad_width(test_args: Args) -> None:
    bad_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            test_args.strip_thickness, test_args.thickness_of_substrate, bad_width)
    with raises(TypeError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            test_args.strip_thickness, test_args.thickness_of_substrate, 100)
