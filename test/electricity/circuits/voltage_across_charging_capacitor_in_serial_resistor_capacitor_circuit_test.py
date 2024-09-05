from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity.circuits import voltage_across_charging_capacitor_in_serial_resistor_capacitor_circuit as rc_node

# Description
## Assert we have 3 Volts applied to 2-Ohm resistor in series with 2 Farads capacitor.
## After 1 Tau seconds capacitor voltage should be 63% of initial voltage. Tau = R * C = 4.

Args = namedtuple("Args", ["V0", "R", "C", "T"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    V0 = Quantity(3 * units.volt)
    R = Quantity(2 * units.ohm)
    C = Quantity(2 * units.farad)
    T = Quantity(4 * units.second)
    return Args(V0=V0, R=R, C=C, T=T)


def test_basic_voltage(test_args: Args) -> None:
    result = rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, test_args.R,
        test_args.T)
    assert_equal(result, 1.896 * units.volt)


def test_bad_voltage(test_args: Args) -> None:
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(Vb, test_args.C, test_args.R, test_args.T)
    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(100, test_args.C, test_args.R, test_args.T)


def test_bad_capacity(test_args: Args) -> None:
    Cb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(test_args.V0, Cb, test_args.R, test_args.T)
    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(test_args.V0, 100, test_args.R, test_args.T)


def test_bad_resistance(test_args: Args) -> None:
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, Rb, test_args.T)
    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, 100, test_args.T)


def test_bad_time(test_args: Args) -> None:
    Tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, test_args.R, Tb)
    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, test_args.R, 100)
