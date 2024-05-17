from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.transmission_lines.coplanar_lines import wave_resistance_of_coplanar_line_first as resistance_law

## The width of the central electrode of the coplanar line is 5 millimeter, and the distance between electrodes is 10 millimeter.
## The effective permittivity is 4. Then the wave resistance is 60.28 ohm.

Args = namedtuple("Args",
    ["effective_permittivity", "distance_between_electrodes", "central_electrode_width"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    effective_permittivity = 4
    distance_between_electrodes = Quantity(10 * units.millimeter)
    central_electrode_width = Quantity(5 * units.millimeter)
    return Args(effective_permittivity=effective_permittivity,
        distance_between_electrodes=distance_between_electrodes,
        central_electrode_width=central_electrode_width)


def test_basic_wave_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_wave_resistance(test_args.effective_permittivity,
        test_args.distance_between_electrodes, test_args.central_electrode_width)
    assert_equal(result, 60.28 * units.ohm)


def test_bad_permittivity(test_args: Args) -> None:
    bad_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance(bad_permittivity,
            test_args.distance_between_electrodes, test_args.central_electrode_width)


def test_bad_distance_between_electrodes(test_args: Args) -> None:
    bad_distance_between_electrodes = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance(test_args.effective_permittivity,
            bad_distance_between_electrodes, test_args.central_electrode_width)
    with raises(TypeError):
        resistance_law.calculate_wave_resistance(test_args.effective_permittivity, 100,
            test_args.central_electrode_width)


def test_bad_central_electrode_width(test_args: Args) -> None:
    bad_central_electrode_width = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance(test_args.effective_permittivity,
            test_args.distance_between_electrodes, bad_central_electrode_width)
    with raises(TypeError):
        resistance_law.calculate_wave_resistance(test_args.effective_permittivity,
            test_args.distance_between_electrodes, 100)
    with raises(ValueError):
        resistance_law.calculate_wave_resistance(test_args.effective_permittivity,
            test_args.central_electrode_width, test_args.central_electrode_width)
