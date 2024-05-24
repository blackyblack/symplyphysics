from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
    prefixes
)

from symplyphysics.laws.electricity import corona_discharge_current_from_voltage as current_section_law

## The gas coefficient is 3.57e-9 [coulomb / (meter^2 * volt)]. The mobility of charged particles is 1.6e-4 [meter^2 / (volt * second)].
## The voltage is 3.5 kilovolt. The discharge voltage is 1.5 kilovolt. Then the current is 1.5 kilovolt.

Args = namedtuple("Args", [
    "gas_coefficient", "mobility_of_charged_particles", "voltage", "corona_discharge_occurrence_voltage",
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    gas_coefficient = Quantity(3.57e-9 * units.coulomb / (units.meter**2 * units.volt))
    mobility_of_charged_particles = Quantity(1.6e-4 * units.meter**2 / (units.volt * units.second))
    voltage = Quantity(3.5 * prefixes.kilo * units.volt)
    corona_discharge_occurrence_voltage = Quantity(1.5 * prefixes.kilo * units.volt)

    return Args(gas_coefficient=gas_coefficient,
        mobility_of_charged_particles=mobility_of_charged_particles,
        voltage=voltage,
        corona_discharge_occurrence_voltage=corona_discharge_occurrence_voltage)


def test_basic_current(test_args: Args) -> None:
    result = current_section_law.calculate_current(
        test_args.gas_coefficient, test_args.mobility_of_charged_particles, test_args.voltage,
        test_args.corona_discharge_occurrence_voltage)
    assert_equal(result, 4 * prefixes.micro * units.ampere)


def test_bad_gas_coefficient(test_args: Args) -> None:
    gas_coefficient = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_section_law.calculate_current(gas_coefficient,
            test_args.mobility_of_charged_particles, test_args.voltage,
            test_args.corona_discharge_occurrence_voltage)
    with raises(TypeError):
        current_section_law.calculate_current(100,
            test_args.mobility_of_charged_particles, test_args.voltage,
            test_args.corona_discharge_occurrence_voltage)


def test_bad_mobility_of_charged_particles(test_args: Args) -> None:
    mobility_of_charged_particles = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_section_law.calculate_current(test_args.gas_coefficient,
            mobility_of_charged_particles, test_args.voltage, test_args.corona_discharge_occurrence_voltage)
    with raises(TypeError):
        current_section_law.calculate_current(test_args.gas_coefficient,
            100, test_args.voltage, test_args.corona_discharge_occurrence_voltage)


def test_bad_voltage(test_args: Args) -> None:
    voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_section_law.calculate_current(test_args.gas_coefficient,
            test_args.mobility_of_charged_particles, voltage, test_args.corona_discharge_occurrence_voltage)
    with raises(TypeError):
        current_section_law.calculate_current(test_args.gas_coefficient,
            test_args.mobility_of_charged_particles, 100, test_args.corona_discharge_occurrence_voltage)


def test_bad_corona_discharge_occurrence_voltage(test_args: Args) -> None:
    corona_discharge_occurrence_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_section_law.calculate_current(test_args.gas_coefficient,
            test_args.mobility_of_charged_particles, test_args.voltage, corona_discharge_occurrence_voltage)
    with raises(TypeError):
        current_section_law.calculate_current(test_args.gas_coefficient,
            test_args.mobility_of_charged_particles, test_args.voltage, 100)
