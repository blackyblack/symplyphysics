"""
Admittance matrix equation
==========================

The admittance matrix is one of the ways to describe a microwave device. The
:math:`Y`-parameters of the device act as elements of this matrix. The matrix equation
relates the input and output currents to the input and output voltages.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Two-port_network#Admittance_parameters_(y-parameters)>`__.
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

input_input_admittance = clone_as_symbol(symbols.admittance, display_symbol="Y_ii", display_latex="Y_\\text{ii}")
"""
Ratio of :attr:`~input_current` to :attr:`~input_voltage` in case of a short circuit at
the output. See :symbols:`admittance`.
"""

input_output_admittance = clone_as_symbol(symbols.admittance, display_symbol="Y_io", display_latex="Y_\\text{io}")
"""
Ratio of :attr:`~input_current` to :attr:`~output_voltage` in case of a short circuit at
the output. See :symbols:`admittance`.
"""

output_input_admittance = clone_as_symbol(symbols.admittance, display_symbol="Y_oi", display_latex="Y_\\text{oi}")
"""
Ratio of :attr:`~output_current` to :attr:`~input_voltage` in case of a short circuit at
the output. See :symbols:`admittance`.
"""

output_output_admittance = clone_as_symbol(symbols.admittance, display_symbol="Y_oo", display_latex="Y_\\text{oo}")
"""
Ratio of :attr:`~output_current` to :attr:`~output_voltage` in case of a short circuit
at the output. See :symbols:`admittance`.
"""

law = Eq(
    Matrix([input_current, output_current]),
    Matrix([[input_input_admittance, input_output_admittance], [output_input_admittance, output_output_admittance]])
    * Matrix([input_voltage, output_voltage]))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(input_current_=input_current, output_current_=output_current)
@validate_output(units.voltage)
def calculate_voltages(
    input_current_: Quantity, output_current_: Quantity, conductivities_: tuple[tuple[Quantity,
    Quantity], tuple[Quantity, Quantity]]
) -> tuple[Quantity, Quantity]:
    for i, cc in enumerate(conductivities_):
        for j, conductivity in enumerate(cc):
            assert_equivalent_dimension(conductivity, f"conductivities_[{i}][{j}]",
                "calculate_voltages", units.conductance)
    result_voltages = solve(law, [input_voltage, output_voltage], dict=True)[0]
    result_input_voltage = result_voltages[input_voltage]
    result_output_voltage = result_voltages[output_voltage]
    substitutions = {
        input_current: input_current_,
        output_current: output_current_,
        input_input_admittance: conductivities_[0][0],
        input_output_admittance: conductivities_[0][1],
        output_input_admittance: conductivities_[1][0],
        output_output_admittance: conductivities_[1][1],
    }
    result_input_voltage = result_input_voltage.subs(substitutions)
    result_output_voltage = result_output_voltage.subs(substitutions)
    return (Quantity(result_input_voltage), Quantity(result_output_voltage))
