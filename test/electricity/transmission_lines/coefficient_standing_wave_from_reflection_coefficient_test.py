from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.transmission_lines import coefficient_standing_wave_from_reflection_coefficient as coefficient_law

# Description
## The reflection coefficient module is 0.2.
## Then the standing wave coefficient is 1.5.

Args = namedtuple("Args", ["reflection_coefficient_module"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    reflection_coefficient_module = 0.2

    return Args(reflection_coefficient_module=reflection_coefficient_module)


def test_basic_coefficient_standing_wave(test_args: Args) -> None:
    result = coefficient_law.calculate_coefficient_standing_wave(
        test_args.reflection_coefficient_module)
    assert_equal(result, 1.5)


def test_bad_reflection_coefficient_module(test_args: Args) -> None:
    bad_reflection_coefficient_module = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_coefficient_standing_wave(bad_reflection_coefficient_module)
    with raises(ValueError):
        coefficient_law.calculate_coefficient_standing_wave(
            -test_args.reflection_coefficient_module)
