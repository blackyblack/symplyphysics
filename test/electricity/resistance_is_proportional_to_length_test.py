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
## According to online calculator
## (https://systemlines.ru/tekhnicheskie-i-vspomogatelnye-materialy/kalkulyator-rascheta-soprotivleniya-provodnika/)
## it's resistance should be 0.02625 Ohm.


@fixture(name="test_args")
def test_args_fixture():
    resistivity = Quantity(0.0175 * units.ohm * (units.milli * units.meter)**2 / units.meter)
    wire_length = Quantity(3 * units.meter)
    cross_section = Quantity(2 * (units.milli * units.meter)**2)
    Args = namedtuple("Args", ["resistivity", "wire_length", "cross_section"])
    return Args(resistivity=resistivity, wire_length=wire_length, cross_section=cross_section)


def test_basic_resistance(test_args):
    result = wire_law.calculate_resistance(test_args.resistivity, test_args.wire_length,
        test_args.cross_section)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result_current = convert_to(result, units.ohm).subs(units.ohm, 1).evalf(6)
    assert result_current == approx(0.02625, 0.001)


def test_bad_resistivity(test_args):
    rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wire_law.calculate_resistance(rb, test_args.wire_length, test_args.cross_section)
    with raises(TypeError):
        wire_law.calculate_resistance(100, test_args.wire_length, test_args.cross_section)


def test_bad_length(test_args):
    lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wire_law.calculate_resistance(test_args.resistivity, lb, test_args.cross_section)
    with raises(TypeError):
        wire_law.calculate_resistance(test_args.resistivity, 100, test_args.cross_section)


def test_bad_cross_section(test_args):
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wire_law.calculate_resistance(test_args.resistivity, test_args.wire_length, cb)
    with raises(TypeError):
        wire_law.calculate_resistance(test_args.resistivity, test_args.wire_length, 100)
