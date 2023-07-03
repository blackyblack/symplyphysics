from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import dissipated_heat_power_is_proportional_to_current_square as joule_lenz_law

# Description
## Assert we are having 1Amp flowing through 2-Ohm resistor.
## According to Joule-Lenz law we should have amount of heat dissipated on this resistor equals to 1^2 * 2 = 2 Watts.


@fixture(name="test_args")
def test_args_fixture():
    C = Quantity(1 * units.ampere)
    R = Quantity(2 * units.ohm)
    Args = namedtuple("Args", ["C", "R"])
    return Args(C=C, R=R)


def test_basic_power(test_args):
    result = joule_lenz_law.calculate_heat_power(test_args.C, test_args.R)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.power)
    result_power = convert_to(result, units.watt).evalf(2)
    assert result_power == approx(2, 0.01)


def test_bad_current(test_args):
    Ib = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        joule_lenz_law.calculate_heat_power(Ib, test_args.R)
    with raises(TypeError):
        joule_lenz_law.calculate_heat_power(100, test_args.R)


def test_bad_resistance(test_args):
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        joule_lenz_law.calculate_heat_power(test_args.C, Rb)
    with raises(TypeError):
        joule_lenz_law.calculate_heat_power(test_args.C, 100)
