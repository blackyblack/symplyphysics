"""
Transmission matrix of lossless transmission line
=================================================

Knowing the length of the transmission line, as well as the surge impedance of
the line and the propagation constant of the signal, it is possible to calculate the
parameters :math:`A, B, C, D` of the transmission matrix of a lossless line.

**Notes:**

#. See :ref:`Transmission matrix`.

**Conditions:**

#. The transmission line is lossless.

..
    TODO: find link
"""

from sympy import Eq, solve, Matrix, I, sin, cos, evaluate
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    dimensionless,
    convert_to_float,
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

propagation_constant = SymbolNew("b", 1 / units.length)
"""
The **propagation constant** is the inverse of the signal :symbols:`wavelength`.
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _phase = propagation_constant * length

law = Eq(
    Matrix([[voltage_voltage_parameter, voltage_current_parameter], [current_voltage_parameter, current_current_parameter]]),
    Matrix([[cos(_phase), I * surge_impedance * sin(_phase)], [I * (1 / surge_impedance) * sin(_phase), cos(_phase)]]))
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
    result_a = result[voltage_voltage_parameter]
    result_b = result[voltage_current_parameter]
    result_c = result[current_voltage_parameter]
    result_d = result[current_current_parameter]
    substitutions = {
        surge_impedance: characteristic_resistance_,
        length: line_length_,
        propagation_constant: constant_propagation_,
    }
    result_a = convert_to_float(Quantity(result_a.subs(substitutions)))
    result_b = Quantity(result_b.subs(substitutions))
    result_c = Quantity(result_c.subs(substitutions))
    result_d = convert_to_float(Quantity(result_d.subs(substitutions)))
    assert_equivalent_dimension(result_b, 'result_b', "calculate_transmission_matrix",
        units.impedance)
    assert_equivalent_dimension(result_c, 'result_c', "calculate_transmission_matrix",
        units.conductance)
    return ((result_a, result_b), (result_c, result_d))
