from sympy import Eq, solve
from symplyphysics import (Symbol, validate_input, validate_output, dimensionless, convert_to_float)

# Description
## When calculating the total gain of a transistor amplifier, it is also worth taking into account
## the gain coefficients of the input and output matching circuits.

## Law is: G = G1 * G2 * G3, where
## G - full gain of the transistor amplifier,
## G1 - gain of the input matching circuit,
## G2 - transistor gain,
## G3 - gain of the output matching circuit.

full_gain = Symbol("full_gain", dimensionless)

gain_of_input_matching_circuit = Symbol("gain_of_input_matching_circuit", dimensionless)
transistor_gain = Symbol("transistor_gain", dimensionless)
gain_of_output_matching_circuit = Symbol("gain_of_output_matching_circuit", dimensionless)

law = Eq(full_gain, gain_of_input_matching_circuit * transistor_gain * gain_of_output_matching_circuit)


@validate_input(gain_of_input_matching_circuit_=gain_of_input_matching_circuit,
    transistor_gain_=transistor_gain,
    gain_of_output_matching_circuit_=gain_of_output_matching_circuit)
@validate_output(full_gain)
def calculate_full_gain(gain_of_input_matching_circuit_: float, transistor_gain_: float,
    gain_of_output_matching_circuit_: float) -> float:
    result_expr = solve(law, full_gain, dict=True)[0][full_gain]
    result_expr = result_expr.subs({
        gain_of_input_matching_circuit: gain_of_input_matching_circuit_,
        transistor_gain: transistor_gain_,
        gain_of_output_matching_circuit: gain_of_output_matching_circuit_,
    })
    return convert_to_float(result_expr)
