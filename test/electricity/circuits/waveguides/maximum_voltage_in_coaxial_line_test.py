from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.circuits.waveguides import maximum_voltage_in_coaxial_line as voltage_law

## The outer diameter is 9 millimeters. The inner diameter is 1.35 millimeters.
## The dielectric breakdown intensity is equal to 300 volt per meter.
## Then maximum voltage in the waveguide will be 2.561 volt.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", ["breakdown_intensity", "outer_diameter", "inner_diameter"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    breakdown_intensity = Quantity(300 * units.volt / units.meter)
    outer_diameter = Quantity(9 * units.millimeter)
    inner_diameter = Quantity(1.35 * units.millimeter)
    return Args(breakdown_intensity=breakdown_intensity,
        outer_diameter=outer_diameter,
        inner_diameter=inner_diameter)


def test_basic_maximum_voltage(test_args: Args) -> None:
    result = voltage_law.calculate_maximum_voltage(test_args.breakdown_intensity,
        test_args.outer_diameter, test_args.inner_diameter)
    assert_equal(result, 2.561 * units.volt)


def test_bad_breakdown_intensity(test_args: Args) -> None:
    bad_breakdown_intensity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_maximum_voltage(bad_breakdown_intensity, test_args.outer_diameter,
            test_args.inner_diameter)
    with raises(TypeError):
        voltage_law.calculate_maximum_voltage(100, test_args.outer_diameter,
            test_args.inner_diameter)


def test_bad_diameter(test_args: Args) -> None:
    bad_diameter = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_maximum_voltage(test_args.breakdown_intensity, bad_diameter,
            test_args.inner_diameter)
    with raises(TypeError):
        voltage_law.calculate_maximum_voltage(test_args.breakdown_intensity, 100,
            test_args.inner_diameter)
    with raises(errors.UnitsError):
        voltage_law.calculate_maximum_voltage(test_args.breakdown_intensity,
            test_args.outer_diameter, bad_diameter)
    with raises(TypeError):
        voltage_law.calculate_maximum_voltage(test_args.breakdown_intensity,
            test_args.outer_diameter, 100)
    with raises(ValueError):
        voltage_law.calculate_maximum_voltage(test_args.breakdown_intensity,
            test_args.inner_diameter, test_args.outer_diameter)
