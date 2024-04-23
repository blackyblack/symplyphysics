from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity.circuits import approximation_of_transmission_coefficient_of_low_pass_filter as coefficient_law

# Description
## The cutoff frequency is 1000 Hertz. The selected frequency is 500 Hertz. The distortion in the bandwidth is 0.5, and the Butterworth
## filter function of order 1 will be used as an approximating function. Then the value of the transfer coefficient will be equal to 0.94.

Args = namedtuple("Args", [
    "filter_function", "frequency", "bandwidth_distortion"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    frequency = Quantity(500 * units.hertz)
    cutoff_frequency = Quantity(1000 * units.hertz)
    order = 1
    filter_function = (frequency / cutoff_frequency)**order
    bandwidth_distortion = 0.5

    return Args(filter_function=filter_function,
        frequency=frequency,
        bandwidth_distortion=bandwidth_distortion)


def test_basic_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_coefficient(test_args.filter_function, test_args.bandwidth_distortion)
    assert_equal(result, 0.94, tolerance=0.01)


def test_bad_bandwidth_distortion(test_args: Args) -> None:
    bandwidth_distortion = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_coefficient(test_args.filter_function, bandwidth_distortion)
