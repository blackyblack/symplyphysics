from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (errors, units, Quantity)
from symplyphysics.laws.thermodynamics import efficiency_factor as efficiency_law

# Description
## Suppose we have a heat engine that receives 2 Joules of heat from the heater, and gives 1.5 Joules to the refrigerator. Then its efficiency should be (2 - 1.5) / 2 = 0.25


@fixture(name="test_args")
def test_args_fixture():
    q_h = Quantity(2 * units.joule)
    q_r = Quantity(1.5 * units.joule)
    Args = namedtuple("Args", ["q_h", "q_r"])
    return Args(q_h=q_h, q_r=q_r)


def test_basic_efficiency(test_args):
    result = efficiency_law.calculate_efficiency_factor(test_args.q_h, test_args.q_r)
    assert result == approx(0.25, 0.001)


def test_bad_efficiency(test_args):
    q_rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        efficiency_law.calculate_efficiency_factor(test_args.q_h, q_rb)
    with raises(errors.UnitsError):
        efficiency_law.calculate_efficiency_factor(q_rb, test_args.q_r)
    with raises(TypeError):
        efficiency_law.calculate_efficiency_factor(test_args.q_h, 1.5)
    with raises(TypeError):
        efficiency_law.calculate_efficiency_factor(2, test_args.q_h)
