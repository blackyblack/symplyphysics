from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity import electric_intensity_in_gas_gap_between_two_electrodes as intensity_law

# Description
## The coordinate of the point between the electrodes, located on the axis from the cathode to the anode, is 1 centimeter.
## The distance between the electrodes is 4 centimeter. The voltage between the electrodes is 350 volt.
## Then the electric intensity at the point between the electrodes is 6.56 kilovolt per meter.

Args = namedtuple("Args",
    ["coordinate", "distance_between_electrodes", "voltage_between_electrodes"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    coordinate = Quantity(1 * units.centimeter)
    distance_between_electrodes = Quantity(4 * units.centimeter)
    voltage_between_electrodes = Quantity(350 * units.volt)

    return Args(coordinate=coordinate,
        distance_between_electrodes=distance_between_electrodes,
        voltage_between_electrodes=voltage_between_electrodes)


def test_basic_electric_intensity(test_args: Args) -> None:
    result = intensity_law.calculate_electric_intensity(test_args.coordinate,
        test_args.distance_between_electrodes, test_args.voltage_between_electrodes)
    assert_equal(result, 6.56 * (prefixes.kilo * units.volt / units.meter))


def test_bad_coordinate(test_args: Args) -> None:
    coordinate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_electric_intensity(coordinate,
            test_args.distance_between_electrodes, test_args.voltage_between_electrodes)
    with raises(TypeError):
        intensity_law.calculate_electric_intensity(100, test_args.distance_between_electrodes,
            test_args.voltage_between_electrodes)
    coordinate = Quantity(5 * units.centimeter)
    with raises(ValueError):
        intensity_law.calculate_electric_intensity(coordinate,
            test_args.distance_between_electrodes, test_args.voltage_between_electrodes)


def test_bad_distance_between_electrodes(test_args: Args) -> None:
    distance_between_electrodes = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_electric_intensity(test_args.coordinate,
            distance_between_electrodes, test_args.voltage_between_electrodes)
    with raises(TypeError):
        intensity_law.calculate_electric_intensity(test_args.coordinate, 100,
            test_args.voltage_between_electrodes)


def test_bad_voltage_between_electrodes(test_args: Args) -> None:
    voltage_between_electrodes = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_electric_intensity(test_args.coordinate,
            test_args.distance_between_electrodes, voltage_between_electrodes)
    with raises(TypeError):
        intensity_law.calculate_electric_intensity(test_args.coordinate,
            test_args.distance_between_electrodes, 100)
