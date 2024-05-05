from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.transmission_lines.coplanar_lines import effective_permettivity_for_coplanar_transmission_line_for_distance_greater_thickness as permittivity_law

## The relative permittivity of the coplanar line dielectric is 4. The width of the central electrode of the coplanar line is 5 millimeter, and
## the thickness of the substrate is 2 millimeters, the distance between electrodes is 10 millimeter.
## Then the effective permittivity is 1.895.

Args = namedtuple("Args", [
    "relative_permittivity", "distance_between_electrodes", "thickness_of_substrate", "central_electrode_width"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 4
    distance_between_electrodes = Quantity(10 * units.millimeter)
    thickness_of_substrate = Quantity(2 * units.millimeter)
    central_electrode_width = Quantity(5 * units.millimeter)
    return Args(relative_permittivity=relative_permittivity,
        distance_between_electrodes=distance_between_electrodes,
        thickness_of_substrate=thickness_of_substrate,
        central_electrode_width=central_electrode_width)


def test_basic_effective_permittivity(test_args: Args) -> None:
    result = permittivity_law.calculate_effective_permittivity(
        test_args.relative_permittivity, test_args.distance_between_electrodes, test_args.thickness_of_substrate,
        test_args.central_electrode_width)
    assert_equal(result, 1.895)


def test_bad_permittivity(test_args: Args) -> None:
    bad_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(bad_permittivity,
            test_args.distance_between_electrodes, test_args.thickness_of_substrate,
            test_args.central_electrode_width)


def test_bad_distance_between_electrodes(test_args: Args) -> None:
    bad_distance_between_electrodes = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            bad_distance_between_electrodes, test_args.thickness_of_substrate, test_args.central_electrode_width)
    with raises(TypeError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            100, test_args.thickness_of_substrate, test_args.central_electrode_width)


def test_bad_thickness_of_substrate(test_args: Args) -> None:
    bad_thickness_of_substrate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            test_args.distance_between_electrodes, bad_thickness_of_substrate, test_args.central_electrode_width)
    with raises(TypeError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            test_args.distance_between_electrodes, 100, test_args.central_electrode_width)
    with raises(ValueError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity, test_args.thickness_of_substrate,
                                                          test_args.distance_between_electrodes, test_args.central_electrode_width)


def test_bad_central_electrode_width(test_args: Args) -> None:
    bad_central_electrode_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            test_args.distance_between_electrodes, test_args.thickness_of_substrate, bad_central_electrode_width)
    with raises(TypeError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            test_args.distance_between_electrodes, test_args.thickness_of_substrate, 100)
    with raises(ValueError):
        permittivity_law.calculate_effective_permittivity(
                test_args.relative_permittivity, test_args.central_electrode_width, test_args.thickness_of_substrate,
                test_args.central_electrode_width)
