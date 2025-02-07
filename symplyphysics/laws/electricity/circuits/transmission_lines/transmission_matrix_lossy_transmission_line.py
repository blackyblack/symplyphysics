"""
Transmission matrix of lossy transmission line
==============================================

The transmission parameters matrix is one of the ways to describe a microwave device.
The :math:`ABCD`-parameters of the device act as elements of this matrix. The matrix
equation relates the input voltage and input current to the output voltage and output
current. Knowing the length and the loss factor of the transmission line, as well as the
characteristic resistance of the line and the constant propagation of signal, it is
possible to calculate the parameters :math:`A, B, C, D` of the transmission matrix of
this line.

**Notes:**

#. See :ref:`Transmission matrix`.

..
    TODO: find link
"""

from sympy import Eq, solve, Matrix, S, I, sinh, cosh
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    dimensionless,
    convert_to,
    symbols,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

voltage_voltage_parameter = SymbolNew("A", dimensionless)
"""
Ratio of input :symbols:`voltage` to output :symbols:`voltage` at idle at the output.
"""

voltage_current_parameter = SymbolNew("B", units.impedance)
"""
Ratio of input :symbols:`voltage` to output :symbols:`current` in case of a short circuit at the output.
"""

current_voltage_parameter = SymbolNew("C", units.conductance)
"""
Ratio of input :symbols:`current` to output :symbols:`voltage` at idle at the output.
"""

current_current_parameter = SymbolNew("D", dimensionless)
"""
Ratio of input :symbols:`current` to output :symbols:`current` in case of a short circuit at the output.
"""

surge_impedance = symbols.surge_impedance
"""
:symbols:`surge_impedance` of the transmission line.
"""

length = symbols.length
"""
:symbols:`length` of the transmission line.
"""

propagation_constant = SymbolNew("b", 1 / units.length)
"""
The **propagation constant** is the inverse of the signal :symbols:`wavelength`.
"""

loss_factor = SymbolNew("a", 1 / units.length)
"""
The **loss factor** shows how many times the transmitted signal weakens per unit
:symbols:`length` of the transmission line.
"""

expression = (loss_factor + I * propagation_constant) * length

law = Eq(
    Matrix([[voltage_voltage_parameter, voltage_current_parameter], [current_voltage_parameter, current_current_parameter]]),
    Matrix([
        [cosh(expression), surge_impedance * sinh(expression)],
        [(1 / surge_impedance) * sinh(expression), cosh(expression)],
    ]),
)
"""
..
    NOTE: Code printing in upcoming PR.
"""


@validate_input(characteristic_resistance_=surge_impedance,
    line_length_=length,
    constant_propagation_=propagation_constant,
    loss_factor_=loss_factor)
def calculate_transmission_matrix(
        characteristic_resistance_: Quantity, line_length_: Quantity,
        constant_propagation_: Quantity,
        loss_factor_: Quantity) -> tuple[tuple[float, Quantity], tuple[Quantity, float]]:
    result = solve(law, [
        voltage_voltage_parameter, voltage_current_parameter, current_voltage_parameter,
        current_current_parameter
    ],
        dict=True)[0]
    result_a = result[voltage_voltage_parameter]
    result_b = result[voltage_current_parameter]
    result_c = result[current_voltage_parameter]
    result_d = result[current_current_parameter]
    substitutions = {
        surge_impedance: characteristic_resistance_,
        length: line_length_,
        propagation_constant: constant_propagation_,
        loss_factor: loss_factor_,
    }
    result_a = convert_to(result_a.subs(substitutions), S.One)
    result_b = Quantity(result_b.subs(substitutions))
    result_c = Quantity(result_c.subs(substitutions))
    result_d = convert_to(result_d.subs(substitutions), S.One)
    assert_equivalent_dimension(result_b, 'result_b', "calculate_transmission_matrix",
        units.impedance)
    assert_equivalent_dimension(result_c, 'result_c', "calculate_transmission_matrix",
        units.conductance)
    return ((result_a, result_b), (result_c, result_d))
