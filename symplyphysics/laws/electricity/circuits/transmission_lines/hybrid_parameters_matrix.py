"""
Hybrid parameters matrix
========================

The **hybrid parameters matrix** is one of the ways to describe a microwave device. The
:math:`H`-parameters of the device act as elements of this matrix. The matrix equation
relates the input voltage and output current to the input current and output voltage.

..
    TODO: find link
"""

from sympy import Eq, solve, Matrix
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    dimensionless,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

input_voltage = clone_as_symbol(symbols.voltage, display_symbol="V_i", display_latex="V_\\text{i}")
"""
Input :symbols:`voltage`.
"""

output_voltage = clone_as_symbol(symbols.voltage, display_symbol="V_o", display_latex="V_\\text{o}")
"""
Output :symbols:`voltage`.
"""

input_current = clone_as_symbol(symbols.current, display_symbol="I_i", display_latex="I_\\text{i}")
"""
Input :symbols:`current`.
"""

output_current = clone_as_symbol(symbols.current, display_symbol="I_o", display_latex="I_\\text{o}")
"""
Output :symbols:`current`.
"""

input_input_parameter = SymbolNew("H_ii", units.impedance, display_latex="H_\\text{ii}")
"""
Ratio of input :symbols:`voltage` to input :symbols:`current` in case of a short circuit
at the output.
"""

input_output_parameter = SymbolNew("H_io", dimensionless, display_latex="H_\\text{io}")
"""
Ratio of input :symbols:`voltage` to output :symbols:`voltage` at idle at the input.
"""

output_input_parameter = SymbolNew("H_oi", dimensionless, display_latex="H_\\text{oi}")
"""
Ratio of output :symbols:`current` to input :symbols:`current` in case of a short circuit
at the output.
"""

output_output_parameter = SymbolNew("H_oo", units.conductance, display_latex="H_\\text{oo}")
"""
Ratio of output :symbols:`current` to output :symbols:`voltage` at idle at the input.
"""

law = Eq(
    Matrix([input_voltage, output_current]),
    Matrix([[input_input_parameter, input_output_parameter], [output_input_parameter, output_output_parameter]])
    * Matrix([input_current, output_voltage]))
"""
..
    NOTE: Code printing in upcoming PR.
"""


@validate_input(
    input_voltage_=input_voltage,
    output_current_=output_current,
)
def calculate_current_and_voltage(
    input_voltage_: Quantity, output_current_: Quantity, parameters_: tuple[tuple[Quantity, float],
    tuple[float, Quantity]]) -> tuple[Quantity, Quantity]:
    assert_equivalent_dimension(parameters_[0][0], "parameters_[0][0]",
        "calculate_current_and_voltage", units.impedance)
    assert_equivalent_dimension(parameters_[1][1], "parameters_[1][1]",
        "calculate_current_and_voltage", units.conductance)
    result = solve(law, [input_current, output_voltage], dict=True)[0]
    result_input_current = result[input_current]
    result_output_voltage = result[output_voltage]
    substitutions = {
        input_voltage: input_voltage_,
        output_current: output_current_,
        input_input_parameter: parameters_[0][0],
        input_output_parameter: parameters_[0][1],
        output_input_parameter: parameters_[1][0],
        output_output_parameter: parameters_[1][1],
    }
    result_input_current = Quantity(result_input_current.subs(substitutions))
    result_output_voltage = Quantity(result_output_voltage.subs(substitutions))
    assert_equivalent_dimension(result_input_current, 'result_input_current',
        "calculate_current_and_voltage", units.current)
    assert_equivalent_dimension(result_output_voltage, 'result_output_voltage',
        "calculate_current_and_voltage", units.voltage)
    return (result_input_current, result_output_voltage)
