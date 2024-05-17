from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.transmission_lines.microstrip_lines import effective_permittivity_for_microstrip_line_from_frequency as permittivity_law

## The relative permittivity of the microstrip line dielectric is 4. The width of the microstrip line is 5 millimeter, and
## the thickness of the substrate is 4 millimeters, the frequency is 2 gigahertz. The effective permittivity without taking into account
## the dependence on frequency is 2.95. Then the effective permittivity is 3.06.

Args = namedtuple("Args", [
    "relative_permittivity", "frequency", "thickness_of_substrate", "width",
    "effective_permittivity_without_frequency"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 4
    frequency = Quantity(2 * prefixes.giga * units.hertz)
    thickness_of_substrate = Quantity(4 * units.millimeter)
    width = Quantity(5 * units.millimeter)
    effective_permittivity_without_frequency = 2.95
    return Args(relative_permittivity=relative_permittivity,
        frequency=frequency,
        thickness_of_substrate=thickness_of_substrate,
        width=width,
        effective_permittivity_without_frequency=effective_permittivity_without_frequency)


def test_basic_effective_permittivity(test_args: Args) -> None:
    result = permittivity_law.calculate_effective_permittivity(
        test_args.relative_permittivity, test_args.frequency, test_args.thickness_of_substrate,
        test_args.width, test_args.effective_permittivity_without_frequency)
    assert_equal(result, 3.06)


def test_bad_permittivities(test_args: Args) -> None:
    bad_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(
            bad_permittivity, test_args.frequency, test_args.thickness_of_substrate,
            test_args.width, test_args.effective_permittivity_without_frequency)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity,
            test_args.frequency, test_args.thickness_of_substrate, test_args.width,
            bad_permittivity)


def test_bad_frequency(test_args: Args) -> None:
    bad_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(
            test_args.relative_permittivity, bad_frequency, test_args.thickness_of_substrate,
            test_args.width, test_args.effective_permittivity_without_frequency)
    with raises(TypeError):
        permittivity_law.calculate_effective_permittivity(
            test_args.relative_permittivity, 100, test_args.thickness_of_substrate, test_args.width,
            test_args.effective_permittivity_without_frequency)


def test_bad_thickness_of_substrate(test_args: Args) -> None:
    bad_thickness_of_substrate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(
            test_args.relative_permittivity, test_args.frequency, bad_thickness_of_substrate,
            test_args.width, test_args.effective_permittivity_without_frequency)
    with raises(TypeError):
        permittivity_law.calculate_effective_permittivity(
            test_args.relative_permittivity, test_args.frequency, 100, test_args.width,
            test_args.effective_permittivity_without_frequency)


def test_bad_width(test_args: Args) -> None:
    bad_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(
            test_args.relative_permittivity, test_args.frequency, test_args.thickness_of_substrate,
            bad_width, test_args.effective_permittivity_without_frequency)
    with raises(TypeError):
        permittivity_law.calculate_effective_permittivity(
            test_args.relative_permittivity, test_args.frequency, test_args.thickness_of_substrate,
            100, test_args.effective_permittivity_without_frequency)
