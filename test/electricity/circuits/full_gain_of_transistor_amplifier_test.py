from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import full_gain_of_transistor_amplifier as gain_law

# Description
## The gain of the input matching circuit is 0.8. The transistor gain is 3. 
## The gain of the output matching circuit is 0.6. 
## Then the full gain of the transistor amplifier is 1.44.

Args = namedtuple("Args", ["gain_of_input_matching_circuit", "transistor_gain", "gain_of_output_matching_circuit"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    gain_of_input_matching_circuit = 0.8
    transistor_gain = 3
    gain_of_output_matching_circuit = 0.6

    return Args(gain_of_input_matching_circuit=gain_of_input_matching_circuit,
        transistor_gain=transistor_gain,
        gain_of_output_matching_circuit=gain_of_output_matching_circuit)


def test_basic_full_gain(test_args: Args) -> None:
    result = gain_law.calculate_full_gain(test_args.gain_of_input_matching_circuit, test_args.transistor_gain,
        test_args.gain_of_output_matching_circuit)
    assert_equal(result, 1.44)


def test_bad_gain_of_input_matching_circuit(test_args: Args) -> None:
    gain_of_input_matching_circuit = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gain_law.calculate_full_gain(gain_of_input_matching_circuit, test_args.transistor_gain,
            test_args.gain_of_output_matching_circuit)


def test_bad_transistor_gain(test_args: Args) -> None:
    transistor_gain = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gain_law.calculate_full_gain(test_args.gain_of_input_matching_circuit, transistor_gain,
            test_args.gain_of_output_matching_circuit)


def test_bad_gain_of_output_matching_circuit(test_args: Args) -> None:
    gain_of_output_matching_circuit = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gain_law.calculate_full_gain(test_args.gain_of_input_matching_circuit, test_args.transistor_gain,
            gain_of_output_matching_circuit)
