from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.circuits.diodes import limit_operating_frequency_of_vacuum_diode as frequency_law

# Description
## The distance between the electrodes of the vacuum diode is 1 centimeter. The voltage between cathode and anode is 17 volt.
## Then the limit operating frequency is 40.76 megahertz.

Args = namedtuple("Args", ["distance_between_electrodes", "voltage"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    distance_between_electrodes = Quantity(1 * units.centimeter)
    voltage = Quantity(17 * units.volt)

    return Args(distance_between_electrodes=distance_between_electrodes, voltage=voltage)


def test_basic_limit_operating_frequency(test_args: Args) -> None:
    result = frequency_law.calculate_limit_operating_frequency(
        test_args.distance_between_electrodes, test_args.voltage)
    assert_equal(result, 40.76 * prefixes.mega * units.hertz)


def test_bad_distance_between_electrodes(test_args: Args) -> None:
    distance_between_electrodes = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_limit_operating_frequency(distance_between_electrodes,
            test_args.voltage)
    with raises(TypeError):
        frequency_law.calculate_limit_operating_frequency(100, test_args.voltage)


def test_bad_voltage(test_args: Args) -> None:
    voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_limit_operating_frequency(test_args.distance_between_electrodes,
            voltage)
    with raises(TypeError):
        frequency_law.calculate_limit_operating_frequency(test_args.distance_between_electrodes,
            100)
