from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import resistance_is_proportional_to_length as wire_law

# Description
## Assert we have 3 meters of copper wire with 2 mm^2 section. Resistivity of copper is 1.75e-8 Ohm*m.
## According to online calculator (https://systemlines.ru/tekhnicheskie-i-vspomogatelnye-materialy/kalkulyator-rascheta-soprotivleniya-provodnika/) it's resistance should be 0,02625 Ohm.

@fixture
def test_args():
    resistivity = Quantity(0.0172 * units.ohm * (units.milli * units.meter)^2 / units.meter)
    wire_length = Quantity(3 * units.meter)
    cross_section = Quantity(2 * (units.milli * units.meter)^2)
    Args = namedtuple("Args", ["Resistivity", "Wire_length", "Cross_section"])
    return Args(resistivity=resistivity, wire_length=wire_length, cross_section=cross_section)


def test_basic_resistance(test_args):
    result = wire_law.calculate_resistance(test_args.Resistivity, test_args.Length, test_args.Section)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.ohm)
    result_current = convert_to(result, units.ohm).subs(units.ohm, 1).evalf(2)
    assert result_current == approx(0.02625, 0.001)

'''
def test_bad_voltage(test_args):
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        ohms_law.calculate_current(Vb, test_args.Resistance)
    with raises(TypeError):
        ohms_law.calculate_current(100, test_args.Resistance)

'''
