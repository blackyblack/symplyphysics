from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.circuits import conductivity_for_rectangular_loop_coupler as conductivity_law

## The standard transmission line conductivity is 20 millisiemens, the power ratio at the output ports equal to 2. 
## Then the values of Y1, Y2, Y3, Y4 impedances are equal, respectively: 14.14 millisiemens, 24.49 millisiemens, 24.49 millisiemens, 14.14 millisiemens.

Args = namedtuple("Args", ["transmission_line_conductivity", "ratio_of_power"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    transmission_line_conductivity = Quantity(20 * prefixes.milli * units.siemens)
    ratio_of_power = 2

    return Args(transmission_line_conductivity=transmission_line_conductivity,
        ratio_of_power=ratio_of_power
        )


def test_basic_conductivities(test_args: Args) -> None:
    result = conductivity_law.calculate_conductivities(test_args.transmission_line_conductivity, test_args.ratio_of_power)
    assert_equal(result[0], 14.14 * prefixes.milli * units.siemens)
    assert_equal(result[1], 24.49 * prefixes.milli * units.siemens)
    assert_equal(result[2], 24.49 * prefixes.milli * units.siemens)
    assert_equal(result[3], 14.14 * prefixes.milli * units.siemens)


def test_bad_transmission_line_conductivity(test_args: Args) -> None:
    bad_transmission_line_conductivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        conductivity_law.calculate_conductivities(bad_transmission_line_conductivity, test_args.ratio_of_power)
    with raises(TypeError):
        conductivity_law.calculate_conductivities(100, test_args.ratio_of_power, test_args)


def test_bad_ratio_of_power(test_args: Args) -> None:
    bad_ratio_of_power = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        conductivity_law.calculate_conductivities(test_args.transmission_line_conductivity, bad_ratio_of_power)
