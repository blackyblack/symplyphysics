from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity, prefixes)
from symplyphysics.laws.electricity import resistance_is_proportional_to_length as wire_law

# Description
## Assert we have 3 meters of copper wire with 2 mm^2 section. Resistivity of copper is 1.75e-8 Ohm*m.
## According to online calculator
## (https://systemlines.ru/tekhnicheskie-i-vspomogatelnye-materialy/kalkulyator-rascheta-soprotivleniya-provodnika/)
## it's resistance should be 0.02625 Ohm.

Args = namedtuple("Args", ["resistivity", "wire_length", "cross_section"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    resistivity = Quantity(0.0175 * units.ohm * (prefixes.milli * units.meter)**2 / units.meter)
    wire_length = Quantity(3 * units.meter)
    cross_section = Quantity(2 * (prefixes.milli * units.meter)**2)
    return Args(resistivity=resistivity, wire_length=wire_length, cross_section=cross_section)


def test_basic_resistance(test_args: Args) -> None:
    result = wire_law.calculate_resistance(test_args.resistivity, test_args.wire_length,
        test_args.cross_section)
    assert_equal(result, 0.02625 * units.ohm)


def test_bad_resistivity(test_args: Args) -> None:
    rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wire_law.calculate_resistance(rb, test_args.wire_length, test_args.cross_section)
    with raises(TypeError):
        wire_law.calculate_resistance(100, test_args.wire_length, test_args.cross_section)


def test_bad_length(test_args: Args) -> None:
    lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wire_law.calculate_resistance(test_args.resistivity, lb, test_args.cross_section)
    with raises(TypeError):
        wire_law.calculate_resistance(test_args.resistivity, 100, test_args.cross_section)


def test_bad_cross_section(test_args: Args) -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wire_law.calculate_resistance(test_args.resistivity, test_args.wire_length, cb)
    with raises(TypeError):
        wire_law.calculate_resistance(test_args.resistivity, test_args.wire_length, 100)
