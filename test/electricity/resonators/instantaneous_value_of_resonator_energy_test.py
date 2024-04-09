from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.resonators import instantaneous_value_of_resonator_energy as energy_law

## The initial energy value in the resonator is 10000 joule. The frequency is 2 * pi * 100e6 radian per second,
## the time is 1 microsecond. The quality factor of the resonator is 10000.
## Then the energy value will be 9391 joules.

Args = namedtuple("Args", ["initial_energy", "time", "angular_frequency", "quality_factor"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    initial_energy = Quantity(10000 * units.joule)
    time = Quantity(1 * prefixes.micro * units.second)
    angular_frequency = Quantity(2 * pi * 100e6 * (units.radian / units.second))
    quality_factor = 10000
    return Args(initial_energy=initial_energy,
        time=time,
        angular_frequency=angular_frequency,
        quality_factor=quality_factor)


def test_basic_instantaneous_energy(test_args: Args) -> None:
    result = energy_law.calculate_instantaneous_energy(test_args.initial_energy, test_args.time,
        test_args.angular_frequency, test_args.quality_factor)
    assert_equal(result, 9391 * units.joule)


def test_bad_initial_energy(test_args: Args) -> None:
    bad_initial_energy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_instantaneous_energy(bad_initial_energy, test_args.time,
            test_args.angular_frequency, test_args.quality_factor)
    with raises(TypeError):
        energy_law.calculate_instantaneous_energy(100, test_args.time, test_args.angular_frequency,
            test_args.quality_factor)


def test_bad_time(test_args: Args) -> None:
    bad_time = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_instantaneous_energy(test_args.initial_energy, bad_time,
            test_args.angular_frequency, test_args.quality_factor)
    with raises(TypeError):
        energy_law.calculate_instantaneous_energy(test_args.initial_energy, 100,
            test_args.angular_frequency, test_args.quality_factor)


def test_bad_angular_frequency(test_args: Args) -> None:
    bad_angular_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_instantaneous_energy(test_args.initial_energy, test_args.time,
            bad_angular_frequency, test_args.quality_factor)
    with raises(TypeError):
        energy_law.calculate_instantaneous_energy(test_args.initial_energy, test_args.time, 100,
            test_args.quality_factor)


def test_bad_quality_factor(test_args: Args) -> None:
    quality_factor = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_instantaneous_energy(test_args.initial_energy, test_args.time,
            test_args.angular_frequency, quality_factor)
