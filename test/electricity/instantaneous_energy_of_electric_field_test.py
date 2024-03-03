from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import instantaneous_energy_of_electric_field as energy_law

# Description
## The inductance in the oscillatory circuit is 1 henry, the maximum current value is 0.02 ampere, and the frequency is 2 * pi * 200 radian per second.
## With an initial phase equal to zero and a time equal to 5e-3 second, the value of the electric field energy will be equal to 2e-4 joules.
## https://remote.misis.ru/courses/168/pages/13-dot-7-primiery-rieshieniia-zadach

Args = namedtuple("Args",
    ["inductance", "maximum_current", "frequency", "time", "initial_phase"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    inductance = Quantity(1 * units.henry)
    maximum_current = Quantity(0.02 * units.ampere)
    frequency = Quantity(2 * pi * 200 * units.radian / units.second)
    time = Quantity(5e-3 * units.second)
    initial_phase = 0
    return Args(inductance=inductance,
        maximum_current=maximum_current,
        frequency=frequency,
        time=time,
        initial_phase=initial_phase)


def test_basic_energy(test_args: Args) -> None:
    result = energy_law.calculate_energy(test_args.inductance, test_args.maximum_current,
        test_args.frequency, test_args.time, test_args.initial_phase)
    assert_equal(result, 2e-4 * units.joule)


def test_bad_inductance(test_args: Args) -> None:
    inductance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_energy(inductance, test_args.maximum_current,
            test_args.frequency, test_args.time, test_args.initial_phase)
    with raises(TypeError):
        energy_law.calculate_energy(100, test_args.maximum_current,
            test_args.frequency, test_args.time, test_args.initial_phase)


def test_bad_maximum_current(test_args: Args) -> None:
    maximum_current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_energy(test_args.inductance, maximum_current,
            test_args.frequency, test_args.time, test_args.initial_phase)
    with raises(TypeError):
        energy_law.calculate_energy(test_args.inductance, 100,
            test_args.frequency, test_args.time, test_args.initial_phase)


def test_bad_frequency(test_args: Args) -> None:
    frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_energy(test_args.inductance, test_args.maximum_current,
            frequency, test_args.time, test_args.initial_phase)
    with raises(TypeError):
        energy_law.calculate_energy(test_args.inductance, test_args.maximum_current, 100,
            test_args.time, test_args.initial_phase)


def test_bad_time(test_args: Args) -> None:
    time = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_energy(test_args.inductance, test_args.maximum_current,
            test_args.frequency, time, test_args.initial_phase)
    with raises(TypeError):
        energy_law.calculate_energy(test_args.inductance, test_args.maximum_current,
            test_args.frequency, 100, test_args.initial_phase)


def test_bad_initial_phase(test_args: Args) -> None:
    initial_phase = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_energy(test_args.inductance, test_args.maximum_current,
            test_args.frequency, test_args.time, initial_phase)
