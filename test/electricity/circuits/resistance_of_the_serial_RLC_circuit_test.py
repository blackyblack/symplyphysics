from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import resistance_of_the_serial_RLC_circuit as resistance_law

# Description
## With a resistor resistance of 3 ohm, an inductive resistance of 7 ohm, a capacitive resistance of 3 ohm,
# the circuit resistance will be 5 ohm.

Args = namedtuple("Args", ["resistance_resistor", "capacitive_resistance", "inductive_resistance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    resistance_resistor = Quantity(3 * units.ohm)
    capacitive_resistance = Quantity(3 * units.ohm)
    inductive_resistance = Quantity(7 * units.ohm)
    return Args(resistance_resistor=resistance_resistor, capacitive_resistance=capacitive_resistance, inductive_resistance=inductive_resistance)


def test_basic_circuit_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_circuit_resistance(test_args.resistance_resistor, test_args.capacitive_resistance,
        test_args.inductive_resistance)
    assert_equal(result, 5 * units.ohm)


def test_bad_resistance_resistor(test_args: Args) -> None:
    resistance_resistor = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_circuit_resistance(resistance_resistor, test_args.capacitive_resistance,
            test_args.inductive_resistance)
    with raises(TypeError):
        resistance_law.calculate_circuit_resistance(100, test_args.capacitive_resistance, test_args.inductive_resistance)


def test_bad_capacitive_resistance(test_args: Args) -> None:
    capacitive_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_circuit_resistance(test_args.resistance_resistor, capacitive_resistance,
            test_args.inductive_resistance)
    with raises(TypeError):
        resistance_law.calculate_circuit_resistance(test_args.resistance_resistor, 100, test_args.inductive_resistance)


def test_bad_inductive_resistance(test_args: Args) -> None:
    inductive_resistance = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        resistance_law.calculate_circuit_resistance(test_args.resistance_resistor, test_args.capacitive_resistance,
            inductive_resistance)
    with raises(TypeError):
        resistance_law.calculate_circuit_resistance(test_args.resistance_resistor, test_args.capacitive_resistance, 100)
