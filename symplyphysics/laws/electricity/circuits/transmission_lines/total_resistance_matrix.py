"""
Impedance matrix equation
=========================

The impedance matrix is one of the ways to describe a microwave device. The
:math:`Z`-parameters of the device act as elements of this matrix. The matrix equation
relates the input and output voltages to the input and output currents.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Two-port_network#Impedance_parameters_(z-parameters)>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
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

input_input_impedance = clone_as_symbol(symbols.electrical_impedance, display_symbol="Z_ii", display_latex="Z_\\text{ii}")
"""
Ratio of :attr:`~input_voltage` to :attr:`~input_current` at idle at the output. See
:symbols:`electrical_impedance`.
"""

input_output_impedance = clone_as_symbol(symbols.electrical_impedance, display_symbol="Z_io", display_latex="Z_\\text{io}")
"""
Ratio of :attr:`~input_voltage` to :attr:`~output_current` at idle at the output. See
:symbols:`electrical_impedance`.
"""

output_input_impedance = clone_as_symbol(symbols.electrical_impedance, display_symbol="Z_oi", display_latex="Z_\\text{oi}")
"""
Ratio of :attr:`~output_voltage` to :attr:`~input_current` at idle at the output. See
:symbols:`electrical_impedance`.
"""

output_output_impedance = clone_as_symbol(symbols.electrical_impedance, display_symbol="Z_oo", display_latex="Z_\\text{oo}")
"""
Ratio of :attr:`~output_voltage` to :attr:`~output_current` at idle at the output. See
:symbols:`electrical_impedance`.
"""

law = Eq(
    Matrix([input_voltage, output_voltage]),
    Matrix([[input_input_impedance, input_output_impedance],
    [output_input_impedance, output_output_impedance]]) * Matrix([input_current, output_current]))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    input_voltage_=input_voltage,
    output_voltage_=output_voltage,
)
@validate_output(units.current)
def calculate_currents(
    input_voltage_: Quantity, output_voltage_: Quantity, impedances_: tuple[tuple[Quantity,
    Quantity], tuple[Quantity, Quantity]]
) -> tuple[Quantity, Quantity]:
    for i, cc in enumerate(impedances_):
        for j, impedance in enumerate(cc):
            assert_equivalent_dimension(impedance, f"impedances_[{i}][{j}]", "calculate_currents",
                units.impedance)
    result_currents = solve(law, [input_current, output_current], dict=True)[0]
    result_input_current = result_currents[input_current]
    result_output_current = result_currents[output_current]
    substitutions = {
        input_voltage: input_voltage_,
        output_voltage: output_voltage_,
        input_input_impedance: impedances_[0][0],
        input_output_impedance: impedances_[0][1],
        output_input_impedance: impedances_[1][0],
        output_output_impedance: impedances_[1][1],
    }
    result_input_current = result_input_current.subs(substitutions)
    result_output_current = result_output_current.subs(substitutions)
    return (Quantity(result_input_current), Quantity(result_output_current))
