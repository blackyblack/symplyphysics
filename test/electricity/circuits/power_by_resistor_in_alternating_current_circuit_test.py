from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import power_by_resistor_in_alternating_current_circuit as power_law

# Description
## The amplitude of the current is 1 ampere, the resistance of the resistor is 100 ohm, the frequency of the alternating current
## is 2 * pi * 1e3 radian per second, the time is 250 microsecond. Then the power is 100 watt.

Args = namedtuple("Args", ["current_amplitude", "resistance", "current_frequency", "time"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    current_amplitude = Quantity(1 * units.ampere)
    resistance = Quantity(100 * units.ohm)
    current_frequency = Quantity(2 * pi * 1e3 * (units.radian / units.second))
    time = Quantity(250 * units.microsecond)
    return Args(current_amplitude=current_amplitude,
        resistance=resistance,
        time=time,
        current_frequency=current_frequency)


def test_basic_power(test_args: Args) -> None:
    result = power_law.calculate_power(test_args.current_amplitude, test_args.resistance,
        test_args.time, test_args.current_frequency)
    assert_equal(result, 100 * units.watt)


def test_bad_current_amplitude(test_args: Args) -> None:
    current_amplitude = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_power(current_amplitude, test_args.resistance, test_args.time,
            test_args.current_frequency)
    with raises(TypeError):
        power_law.calculate_power(100, test_args.resistance, test_args.time,
            test_args.current_frequency)


def test_bad_resistance(test_args: Args) -> None:
    resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_power(test_args.current_amplitude, resistance, test_args.time,
            test_args.current_frequency)
    with raises(TypeError):
        power_law.calculate_power(test_args.current_amplitude, 100, test_args.time,
            test_args.current_frequency)


def test_bad_time(test_args: Args) -> None:
    time = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_power(test_args.current_amplitude, test_args.resistance, time,
            test_args.current_frequency)
    with raises(TypeError):
        power_law.calculate_power(test_args.current_amplitude, test_args.resistance, 100,
            test_args.current_frequency)


def test_bad_current_frequency(test_args: Args) -> None:
    current_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_power(test_args.current_amplitude, test_args.resistance, test_args.time,
            current_frequency)
    with raises(TypeError):
        power_law.calculate_power(test_args.current_amplitude, test_args.resistance, test_args.time,
            100)
