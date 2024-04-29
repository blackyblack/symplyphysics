from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.transmission_lines import ratio_of_average_power_to_incident_power_from_standing_wave_coefficient as coefficient_law

# Description
## The average power delivered to the load is 3 watt. The incident power is 4 watt. Then the standing wave coefficient is 3.

Args = namedtuple("Args", ["incident_power", "average_power"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    incident_power = Quantity(4 * units.watt)
    average_power = Quantity(3 * units.watt)

    return Args(incident_power=incident_power, average_power=average_power)


def test_basic_standing_wave_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_standing_wave_coefficient(test_args.incident_power, test_args.average_power)
    assert_equal(result, 3)


def test_bad_powers(test_args: Args) -> None:
    bad_power = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_standing_wave_coefficient(bad_power, test_args.average_power)
    with raises(TypeError):
        coefficient_law.calculate_standing_wave_coefficient(100, test_args.average_power)
    with raises(errors.UnitsError):
        coefficient_law.calculate_standing_wave_coefficient(test_args.incident_power, bad_power)
    with raises(TypeError):
        coefficient_law.calculate_standing_wave_coefficient(test_args.incident_power, 100)
    with raises(ValueError):
        coefficient_law.calculate_standing_wave_coefficient(test_args.average_power, test_args.incident_power)
