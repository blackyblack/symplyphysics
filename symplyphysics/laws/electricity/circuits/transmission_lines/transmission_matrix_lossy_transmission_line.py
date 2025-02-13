"""
Transmission matrix of lossy transmission line
==============================================

The transmission parameters matrix is one of the ways to describe a microwave device.
The :math:`ABCD`-parameters of the device act as elements of this matrix. The matrix
equation relates the input voltage and input current to the output voltage and output
current. Knowing the length and the loss factor of the transmission line, as well as the
surge impedance of the line and the constant propagation of signal, it is possible to
calculate the parameters :math:`A, B, C, D` of the transmission matrix of this line.

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
Ratio of input :symbols:`voltage` to output :symbols:`current` in case of a short
circuit at the output.
"""

current_voltage_parameter = SymbolNew("C", units.conductance)
"""
Ratio of input :symbols:`current` to output :symbols:`voltage` at idle at the output.
"""

current_current_parameter = SymbolNew("D", dimensionless)
"""
Ratio of input :symbols:`current` to output :symbols:`current` in case of a short
circuit at the output.
"""

surge_impedance = symbols.surge_impedance
"""
:symbols:`surge_impedance` of the transmission line.
"""

length = symbols.length
"""
:symbols:`length` of the transmission line.
"""

propagation_constant = symbols.propagation_constant
"""
:symbols:`propagation_constant`.
"""

_expression = propagation_constant * length

law = Eq(
    Matrix([[voltage_voltage_parameter, voltage_current_parameter], [current_voltage_parameter, current_current_parameter]]),
    Matrix([
        [cosh(_expression), surge_impedance * sinh(_expression)],
        [(1 / surge_impedance) * sinh(_expression), cosh(_expression)],
    ]),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(characteristic_resistance_=surge_impedance,
    line_length_=length,
    constant_propagation_=propagation_constant)
def calculate_transmission_matrix(
        characteristic_resistance_: Quantity, line_length_: Quantity,
        constant_propagation_: Quantity) -> tuple[tuple[float, Quantity], tuple[Quantity, float]]:
    result = solve(law, [
        voltage_voltage_parameter, voltage_current_parameter, current_voltage_parameter,
        current_current_parameter
    ],
        dict=True)[0]
    result_vv = result[voltage_voltage_parameter]
    result_vc = result[voltage_current_parameter]
    result_cv = result[current_voltage_parameter]
    result_cc = result[current_current_parameter]
    substitutions = {
        surge_impedance: characteristic_resistance_,
        length: line_length_,
        propagation_constant: constant_propagation_,
    }
    result_vv = convert_to(result_vv.subs(substitutions), S.One)
    result_vc = Quantity(result_vc.subs(substitutions))
    result_cv = Quantity(result_cv.subs(substitutions))
    result_cc = convert_to(result_cc.subs(substitutions), S.One)
    assert_equivalent_dimension(result_vc, 'result_vc', "calculate_transmission_matrix",
        units.impedance)
    assert_equivalent_dimension(result_cv, 'result_cv', "calculate_transmission_matrix",
        units.conductance)
    return ((result_vv, result_vc), (result_cv, result_cc))
