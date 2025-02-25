"""
Transmission matrix equation
============================

The transmission matrix equation relates the input voltage and input current to the
output voltage and output current.

**Notes:**

#. See :ref:`Transmission Matrix <transmission_matrix_def>`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Two-port_network#ABCD-parameters>`__.

..
    TODO: rename file
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    dimensionless,
    symbols,
    clone_as_symbol,
    Symbol,
    Matrix,
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

voltage_voltage_parameter = Symbol("A", dimensionless)
"""
Ratio of input :symbols:`voltage` to output :symbols:`voltage` at idle at the output.
"""

voltage_current_parameter = Symbol("B", units.impedance)
"""
Ratio of input :symbols:`voltage` to output :symbols:`current` in case of a short
circuit at the output.
"""

current_voltage_parameter = Symbol("C", units.conductance)
"""
Ratio of input :symbols:`current` to output :symbols:`voltage` at idle at the output.
"""

current_current_parameter = Symbol("D", dimensionless)
"""
Ratio of input :symbols:`current` to output :symbols:`current` in case of a short
circuit at the output.
"""

law = Eq(
    Matrix([input_voltage, input_current]),
    Matrix([[voltage_voltage_parameter, voltage_current_parameter], [current_voltage_parameter, current_current_parameter]])
    * Matrix([output_voltage, output_current]))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    input_voltage_=input_voltage,
    input_current_=input_current,
)
def calculate_current_and_voltage(
    input_voltage_: Quantity, input_current_: Quantity, parameters_: tuple[tuple[float, Quantity],
    tuple[Quantity, float]]) -> tuple[Quantity, Quantity]:
    assert_equivalent_dimension(parameters_[0][1], "parameters_[0][1]",
        "calculate_current_and_voltage", units.impedance)
    assert_equivalent_dimension(parameters_[1][0], "parameters_[1][0]",
        "calculate_current_and_voltage", units.conductance)
    result = solve(law, [output_voltage, output_current], dict=True)[0]
    result_output_current = result[output_current]
    result_output_voltage = result[output_voltage]
    substitutions = {
        input_voltage: input_voltage_,
        input_current: input_current_,
        voltage_voltage_parameter: parameters_[0][0],
        voltage_current_parameter: parameters_[0][1],
        current_voltage_parameter: parameters_[1][0],
        current_current_parameter: parameters_[1][1],
    }
    result_output_current = Quantity(result_output_current.subs(substitutions))
    result_output_voltage = Quantity(result_output_voltage.subs(substitutions))
    assert_equivalent_dimension(result_output_current, 'result_output_current',
        "calculate_current_and_voltage", units.current)
    assert_equivalent_dimension(result_output_voltage, 'result_output_voltage',
        "calculate_current_and_voltage", units.voltage)
    return (result_output_voltage, result_output_current)
