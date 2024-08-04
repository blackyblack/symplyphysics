from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits.transmission_lines import reflection_coefficient_from_load_and_transmission_line_impedances as coefficient_law

# Description
## The load impedance is 120 ohm, and the characteristic impedance of the transmission line is 50 ohm.
## Then the reflection coefficient is 0.412.

Args = namedtuple("Args", ["load_impedance", "characteristic_impedance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    load_impedance = Quantity(120 * units.ohm)
    characteristic_impedance = Quantity(50 * units.ohm)

    return Args(load_impedance=load_impedance, characteristic_impedance=characteristic_impedance)


def test_basic_reflection_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_reflection_coefficient(test_args.load_impedance,
        test_args.characteristic_impedance)
    assert_equal(result, 0.412)


def test_bad_powers(test_args: Args) -> None:
    bad_power = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_reflection_coefficient(bad_power,
            test_args.characteristic_impedance)
    with raises(TypeError):
        coefficient_law.calculate_reflection_coefficient(100, test_args.characteristic_impedance)
    with raises(errors.UnitsError):
        coefficient_law.calculate_reflection_coefficient(test_args.load_impedance, bad_power)
    with raises(TypeError):
        coefficient_law.calculate_reflection_coefficient(test_args.load_impedance, 100)
