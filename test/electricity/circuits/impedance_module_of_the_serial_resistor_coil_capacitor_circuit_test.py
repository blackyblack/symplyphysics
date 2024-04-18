from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import impedance_module_of_the_serial_resistor_coil_capacitor_circuit as resistance_law

# Description
## With a resistor resistance of 3 ohm, an inductive reactance of 7 ohm, a capacitive reactance of 3 ohm,
## the circuit impedance module will be 5 ohm.

Args = namedtuple("Args", ["resistance_resistor", "capacitive_reactance", "inductive_reactance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    resistance_resistor = Quantity(3 * units.ohm)
    capacitive_reactance = Quantity(3 * units.ohm)
    inductive_reactance = Quantity(7 * units.ohm)
    return Args(resistance_resistor=resistance_resistor, capacitive_reactance=capacitive_reactance, inductive_reactance=inductive_reactance)


def test_basic_circuit_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_circuit_impedance_module(test_args.resistance_resistor, test_args.capacitive_reactance,
        test_args.inductive_reactance)
    assert_equal(result, 5 * units.ohm)


def test_bad_reactance_resistor(test_args: Args) -> None:
    bad_reactance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_circuit_impedance_module(bad_reactance, test_args.capacitive_reactance,
            test_args.inductive_reactance)
    with raises(TypeError):
        resistance_law.calculate_circuit_impedance_module(100, test_args.capacitive_reactance, test_args.inductive_reactance)
    with raises(errors.UnitsError):
        resistance_law.calculate_circuit_impedance_module(test_args.resistance_resistor, bad_reactance,
            test_args.inductive_reactance)
    with raises(TypeError):
        resistance_law.calculate_circuit_impedance_module(test_args.resistance_resistor, 100, test_args.inductive_reactance)
    with raises(errors.UnitsError):
        resistance_law.calculate_circuit_impedance_module(test_args.resistance_resistor, test_args.capacitive_reactance,
            bad_reactance)
    with raises(TypeError):
        resistance_law.calculate_circuit_impedance_module(test_args.resistance_resistor, test_args.capacitive_reactance, 100)
