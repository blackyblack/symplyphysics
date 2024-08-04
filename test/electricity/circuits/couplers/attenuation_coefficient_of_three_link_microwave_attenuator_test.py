from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits.couplers import attenuation_coefficient_of_three_link_microwave_attenuator as coefficient_law

# Description
## The resistance of the first resistor is 100 ohms, the resistance of the second resistor is 200 ohms.
## Then the attenuation coefficient will be 2.62.

Args = namedtuple("Args", ["first_resistance", "second_resistance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    first_resistance = Quantity(100 * units.ohm)
    second_resistance = Quantity(200 * units.ohm)

    return Args(first_resistance=first_resistance, second_resistance=second_resistance)


def test_basic_attenuation_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_attenuation_coefficient(test_args.first_resistance,
        test_args.second_resistance)
    assert_equal(result, 2.62)


def test_bad_resistance(test_args: Args) -> None:
    bad_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(bad_resistance,
            test_args.second_resistance)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(100, test_args.second_resistance)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.first_resistance,
            bad_resistance)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.first_resistance, 100)
