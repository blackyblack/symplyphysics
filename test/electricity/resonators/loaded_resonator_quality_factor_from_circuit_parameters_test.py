from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.resonators import loaded_resonator_quality_factor_from_circuit_parameters as factor_law

# Description
## The resistance in the oscillatory circuit is 50 ohm, the inductance is 25 nanohenry. The oscillation frequency is 1 kilohertz.
## The resistance in the external circuit is 20 ohm. Then the quality factor will be equal to 45472.

Args = namedtuple("Args", ["resistance", "inductance", "frequency", "load_resistance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    resistance = Quantity(50 * units.ohm)
    inductance = Quantity(25 * prefixes.nano * units.henry)
    frequency = Quantity(2 * prefixes.kilo * units.hertz)
    load_resistance = Quantity(20 * units.ohm)

    return Args(resistance=resistance,
        inductance=inductance,
        frequency=frequency,
        load_resistance=load_resistance)


def test_basic_loaded_resonator_quality_factor(test_args: Args) -> None:
    result = factor_law.calculate_quality_factor(test_args.resistance, test_args.inductance, test_args.frequency, test_args.load_resistance)
    assert_equal(result, 45472)


def test_bad_resistances(test_args: Args) -> None:
    resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(resistance, test_args.inductance, test_args.frequency, test_args.load_resistance)
    with raises(TypeError):
        factor_law.calculate_quality_factor(100, test_args.inductance, test_args.frequency, test_args.load_resistance)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.resistance, test_args.inductance, test_args.frequency, resistance)
    with raises(TypeError):
        factor_law.calculate_quality_factor(test_args.resistance, test_args.inductance, test_args.frequency, 100)


def test_bad_inductance(test_args: Args) -> None:
    inductance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.resistance, inductance, test_args.frequency, test_args.load_resistance)
    with raises(TypeError):
        factor_law.calculate_quality_factor(test_args.resistance, 100, test_args.frequency, test_args.load_resistance)


def test_bad_frequency(test_args: Args) -> None:
    frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.resistance, test_args.inductance, frequency, test_args.load_resistance)
    with raises(TypeError):
            factor_law.calculate_quality_factor(test_args.resistance, test_args.inductance, 100, test_args.load_resistance)
