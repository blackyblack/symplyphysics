"""
Scattering matrix to transmission matrix
========================================

Knowing the parameters of the scattering matrix, it is possible to determine the
parameters of the transmission matrix.

**Notes:**

#. See :ref:`Transmission matrix`.
#. See :ref:`Scattering matrix`.
#. See `this site <https://www.microwaves101.com/encyclopedias/s-parameters>`__ for more
   information about :math:`S`-parameters.
#. See this `Wikipedia page <https://en.wikipedia.org/wiki/Scattering_parameters#Two-port_S-parameters>`__
   for more information about :math:`S`-parameters.

..
    TODO: find link
    TODO: make laws for the S-parameters, see second link above.
"""

from sympy import Eq, solve, Matrix, evaluate
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    dimensionless,
    convert_to_float,
    symbols,
    Symbol
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

voltage_voltage_parameter = Symbol("A", dimensionless)
"""
Ratio of input :symbols:`voltage` to output :symbols:`voltage` at idle at the output.
"""

voltage_current_parameter = Symbol("B", units.impedance)
"""
Ratio of input :symbols:`voltage` to output :symbols:`current` in case of a short circuit at the output.
"""

current_voltage_parameter = Symbol("C", units.conductance)
"""
Ratio of input :symbols:`current` to output :symbols:`voltage` at idle at the output.
"""

current_current_parameter = Symbol("D", dimensionless)
"""
Ratio of input :symbols:`current` to output :symbols:`current` in case of a short circuit at the output.
"""

surge_impedance = symbols.surge_impedance
"""
:symbols:`surge_impedance`.
"""

input_voltage_reflection_coefficient = Symbol("S_ii", dimensionless, display_latex="S_\\text{ii}")
"""
Input port :symbols:`voltage` :symbols:`reflection_coefficient`.
"""

reverse_voltage_gain = Symbol("S_io", dimensionless, display_latex="S_\\text{io}")
"""
Reverse :symbols:`voltage` :symbols:`circuit_gain`.
"""

forward_voltage_gain = Symbol("S_oi", dimensionless, display_latex="S_\\text{oi}")
"""
Forward :symbols:`voltage` :symbols:`circuit_gain`.
"""

output_voltage_reflection_coefficient = Symbol("S_oo", dimensionless, display_latex="S_\\text{oo}")
"""
Output port :symbols:`voltage` :symbols:`reflection_coefficient`.
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _expression_A = ((1 + input_voltage_reflection_coefficient) * (1 - output_voltage_reflection_coefficient) +
        reverse_voltage_gain * forward_voltage_gain) / (2 * forward_voltage_gain)
    _expression_B = surge_impedance * ((1 + input_voltage_reflection_coefficient) *
        (1 + output_voltage_reflection_coefficient) - reverse_voltage_gain * forward_voltage_gain) / (2 *
        forward_voltage_gain)
    _expression_C = (1 / surge_impedance) * ((1 - input_voltage_reflection_coefficient) *
        (1 - output_voltage_reflection_coefficient) - reverse_voltage_gain * forward_voltage_gain) / (2 *
        forward_voltage_gain)
    _expression_D = ((1 - input_voltage_reflection_coefficient) * (1 + output_voltage_reflection_coefficient) +
        reverse_voltage_gain * forward_voltage_gain) / (2 * forward_voltage_gain)

law = Eq(
    Matrix([[voltage_voltage_parameter, voltage_current_parameter], [current_voltage_parameter, current_current_parameter]]),
    Matrix([[_expression_A, _expression_B], [_expression_C, _expression_D]]),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(characteristic_resistance_=surge_impedance,)
def calculate_transmission_matrix(
    characteristic_resistance_: Quantity, parameters_: tuple[tuple[float, float], tuple[float,
    float]]
) -> tuple[tuple[Quantity, float], tuple[float, Quantity]]:
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
        input_voltage_reflection_coefficient: parameters_[0][0],
        reverse_voltage_gain: parameters_[0][1],
        forward_voltage_gain: parameters_[1][0],
        output_voltage_reflection_coefficient: parameters_[1][1],
    }
    result_a = convert_to_float(result_a.subs(substitutions))
    result_b = Quantity(result_b.subs(substitutions))
    result_c = Quantity(result_c.subs(substitutions))
    result_d = convert_to_float(result_d.subs(substitutions))
    assert_equivalent_dimension(result_b, 'result_b', "calculate_transmission_matrix",
        units.impedance)
    assert_equivalent_dimension(result_c, 'result_c', "calculate_transmission_matrix",
        units.conductance)
    return ((result_a, result_b), (result_c, result_d))
