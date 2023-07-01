from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import period_of_spring_from_mass as spring_period

# Description
## According to https://ncalculators.com/mechanical/simple-harmonic-motion-time-period-calculator.htm,
## for object with mass = 10 kg and spring constant = 2.5 N/m, period should be 12.57 seconds.


@fixture(name="test_args")
def test_args_fixture():
    k = Quantity(2.5 * units.newton / units.meter)
    m = Quantity(10 * units.kilogram)
    Args = namedtuple("Args", ["k", "m"])
    return Args(k=k, m=m)


def test_basic_period(test_args):
    result = spring_period.calculate_period(test_args.k, test_args.m)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.time)
    result_period = convert_to(result, units.second).subs(units.second, 1).evalf(2)
    assert result_period == approx(12.57, 0.01)


def test_bad_elascticity(test_args):
    kb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        spring_period.calculate_period(kb, test_args.m)
    with raises(TypeError):
        spring_period.calculate_period(100, test_args.m)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        spring_period.calculate_period(test_args.k, mb)
    with raises(TypeError):
        spring_period.calculate_period(test_args.k, 100)
