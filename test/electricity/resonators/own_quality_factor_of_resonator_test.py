from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.resonators import own_quality_factor_of_resonator as factor_law

# Description
## The resistance in the oscillatory circuit is 50 ohm, the inductance is 25 nanohenry. The oscillation frequency is 2 kilohertz.
## Then the quality factor will be equal to 159154.

Args = namedtuple("Args", ["resistance", "inductance", "frequency"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    resistance = Quantity(50 * units.ohm)
    inductance = Quantity(25 * prefixes.nano * units.henry)
    frequency = Quantity(2 * prefixes.kilo * units.hertz)

    return Args(resistance=resistance, inductance=inductance, frequency=frequency)


def test_basic_quality_factor(test_args: Args) -> None:
    result = factor_law.calculate_quality_factor(test_args.resistance, test_args.inductance,
        test_args.frequency)
    assert_equal(result, 159154)


def test_bad_resistance(test_args: Args) -> None:
    resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(resistance, test_args.inductance, test_args.frequency)
    with raises(TypeError):
        factor_law.calculate_quality_factor(100, test_args.inductance, test_args.frequency)


def test_bad_inductance(test_args: Args) -> None:
    inductance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.resistance, inductance, test_args.frequency)
    with raises(TypeError):
        factor_law.calculate_quality_factor(test_args.resistance, 100, test_args.frequency)


def test_bad_frequency(test_args: Args) -> None:
    frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.resistance, test_args.inductance, frequency)
    with raises(TypeError):
        factor_law.calculate_quality_factor(test_args.resistance, test_args.inductance, 100)
