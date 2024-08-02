from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.laws.thermodynamics import efficiency_of_heat_engine as efficiency_law

# Description
## Suppose we have a heat engine that receives 2 Joules of heat from the heater, and gives 1.5 Joules to the refrigerator. Then its efficiency should be (2 - 1.5) / 2 = 0.25

Args = namedtuple("Args", ["q_h", "q_r"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q_h = Quantity(2 * units.joule)
    q_r = Quantity(1.5 * units.joule)
    return Args(q_h=q_h, q_r=q_r)


def test_basic_efficiency(test_args: Args) -> None:
    result = efficiency_law.calculate_efficiency_factor(test_args.q_h, test_args.q_r)
    assert_equal(result, 0.25)


def test_bad_efficiency(test_args: Args) -> None:
    q_rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        efficiency_law.calculate_efficiency_factor(test_args.q_h, q_rb)
    with raises(errors.UnitsError):
        efficiency_law.calculate_efficiency_factor(q_rb, test_args.q_r)
    with raises(TypeError):
        efficiency_law.calculate_efficiency_factor(test_args.q_h, 1.5)
    with raises(TypeError):
        efficiency_law.calculate_efficiency_factor(2, test_args.q_h)
