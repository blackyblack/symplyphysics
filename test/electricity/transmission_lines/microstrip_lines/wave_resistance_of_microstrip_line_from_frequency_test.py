from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.transmission_lines.microstrip_lines import wave_resistance_of_microstrip_line_from_frequency as resistance_law

## The effective permittivity without frequency of the microstrip line is 2.5, and the effective permittivity of the microstrip line taking into
## account the dependence on frequency is 3. The wave resistance without frequency is 75 Ohm. Then the wave resistance of the microstrip line
## taking into account the dependence on frequency will be equal to 91.29 ohm.

Args = namedtuple("Args", ["wave_resistance_without_frequency", "effective_permittivity", "effective_permittivity_without_frequency"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    wave_resistance_without_frequency = Quantity(75 * units.ohm)
    effective_permittivity = 3
    effective_permittivity_without_frequency = 2.5
    return Args(wave_resistance_without_frequency=wave_resistance_without_frequency,
        effective_permittivity=effective_permittivity,
        effective_permittivity_without_frequency=effective_permittivity_without_frequency)


def test_basic_wave_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_wave_resistance(test_args.wave_resistance_without_frequency,
        test_args.effective_permittivity, test_args.effective_permittivity_without_frequency)
    assert_equal(result, 91.29 * units.ohm)


def test_bad_wave_resistance_without_frequency(test_args: Args) -> None:
    bad_wave_resistance_without_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance(bad_wave_resistance_without_frequency, test_args.effective_permittivity,
            test_args.effective_permittivity_without_frequency)


def test_bad_effective_permittivities(test_args: Args) -> None:
    bad_effective_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance(test_args.wave_resistance_without_frequency, bad_effective_permittivity,
            test_args.effective_permittivity_without_frequency)
    with raises(errors.UnitsError):
        resistance_law.calculate_wave_resistance(test_args.wave_resistance_without_frequency, test_args.effective_permittivity,
            bad_effective_permittivity)
