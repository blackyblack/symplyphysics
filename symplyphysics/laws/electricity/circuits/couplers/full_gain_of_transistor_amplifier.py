"""
Total gain of transistor amplifier
==================================

The total gain of a transistor amplifier depends on the gain of the input and output
matching circuits and the transistor gain.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

total_gain = clone_as_symbol(symbols.circuit_gain)
"""
Total gain of the transistor amplifier. See :symbols:`circuit_gain`.
"""

input_circuit_gain = clone_as_symbol(symbols.circuit_gain, display_symbol="gain_i", display_latex="\\text{gain}_\\text{i}")
"""
Input matching :symbols:`circuit_gain`.
"""

transistor_gain = clone_as_symbol(symbols.circuit_gain, display_symbol="gain_t", display_latex="\\text{gain}_\\text{t}")
"""
Transistor gain. See :symbols:`circuit_gain`.
"""

output_circuit_gain = clone_as_symbol(symbols.circuit_gain, display_symbol="gain_o", display_latex="\\text{gain}_\\text{o}")
"""
Output matching :symbols:`circuit_gain`.
"""

law = Eq(total_gain, input_circuit_gain * transistor_gain * output_circuit_gain)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(gain_of_input_matching_circuit_=input_circuit_gain,
    transistor_gain_=transistor_gain,
    gain_of_output_matching_circuit_=output_circuit_gain)
@validate_output(total_gain)
def calculate_full_gain(gain_of_input_matching_circuit_: float, transistor_gain_: float,
    gain_of_output_matching_circuit_: float) -> float:
    result_expr = solve(law, total_gain, dict=True)[0][total_gain]
    result_expr = result_expr.subs({
        input_circuit_gain: gain_of_input_matching_circuit_,
        transistor_gain: transistor_gain_,
        output_circuit_gain: gain_of_output_matching_circuit_,
    })
    return convert_to_float(result_expr)
