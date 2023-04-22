from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import energy_accumulated_in_capacitor_from_capacitance_and_voltage as capacitor_law

# Description
## Assert we have 0.22uF capacitor charged to 10 volts.
## According to law we should have amount of energy accumulated in this capacitor equals to 0.22 * 0.001 * 10**2 / 2 = 0.011 Joules.


@fixture
def test_args():
    Capacitance = Quantity(0.00022 * units.farad)
    Voltage = Quantity(10 * units.volt)
    Args = namedtuple("Args", ["Capacitance", "Voltage"])
    return Args(Capacitance=Capacitance, Voltage=Voltage)


def test_basic_energy(test_args):
    result = capacitor_law.calculate_accumulated_energy(test_args.Capacitance,
                                                        test_args.Voltage)
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
                                                     units.energy)
    result_power = convert_to(result, units.joule).subs(units.joule, 1).evalf(5)
    assert result_power == approx(0.011, 0.00001)


def test_bad_capacitance(test_args):
    Cb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        capacitor_law.calculate_accumulated_energy(Cb, test_args.Voltage)
    with raises(TypeError):
        capacitor_law.calculate_accumulated_energy(100, test_args.Voltage)


def test_bad_voltage(test_args):
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        capacitor_law.calculate_accumulated_energy(test_args.Capacitance, Vb)
    with raises(TypeError):
        capacitor_law.calculate_accumulated_energy(test_args.Capacitance, 100)
