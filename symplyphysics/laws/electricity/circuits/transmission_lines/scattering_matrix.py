"""
Scattering matrix equation
==========================

The scattering matrix is one of the ways to describe a microwave device. The
:math:`S`-parameters of the device act as elements of this matrix. The matrix equation
relates the waves reflected from the input and output to the waves incident on the input
and output.

**Notes:**

#. For definitions of :math:`a` and :math:`b` presented below, refer
   `here (Wikipedia) <https://en.wikipedia.org/wiki/Scattering_parameters#A_definition>`__.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Scattering_parameters#Two-port_S-parameters>`__.

..
    TODO: rename file
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    dimensionless,
    Symbol,
    Matrix,
)

input_reflected_power_wave = Symbol("b_i", sqrt(units.power), display_latex="b_\\text{i}")
"""
Outgoing ("reflected") power wave at input port.
"""

output_reflected_power_wave = Symbol("b_o", sqrt(units.power), display_latex="b_\\text{o}")
"""
Outgoing ("reflected") power wave at output port.
"""

input_incident_power_wave = Symbol("a_i", sqrt(units.power), display_latex="a_\\text{i}")
"""
Incident power wave at input port.
"""

output_incident_power_wave = Symbol("a_o", sqrt(units.power), display_latex="a_\\text{o}")
"""
Incident power wave at output port.
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

_scattering_matrix = Matrix([
    [input_voltage_reflection_coefficient, reverse_voltage_gain],
    [forward_voltage_gain, output_voltage_reflection_coefficient],
])

law = Eq(
    Matrix([input_reflected_power_wave, output_reflected_power_wave]),
    _scattering_matrix * Matrix([input_incident_power_wave, output_incident_power_wave])
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(input_reflected_wave_=input_reflected_power_wave,
    output_reflected_wave_=output_reflected_power_wave)
@validate_output(sqrt(units.power))
def calculate_waves(
        input_reflected_wave_: Quantity, output_reflected_wave_: Quantity,
        parameters_: tuple[tuple[float, float], tuple[float, float]]) -> tuple[Quantity, Quantity]:
    result = solve(law, [input_incident_power_wave, output_incident_power_wave], dict=True)[0]
    result_input_wave = result[input_incident_power_wave]
    result_output_wave = result[output_incident_power_wave]
    substitutions = {
        input_reflected_power_wave: input_reflected_wave_,
        output_reflected_power_wave: output_reflected_wave_,
        input_voltage_reflection_coefficient: parameters_[0][0],
        reverse_voltage_gain: parameters_[0][1],
        forward_voltage_gain: parameters_[1][0],
        output_voltage_reflection_coefficient: parameters_[1][1],
    }
    result_input_wave = result_input_wave.subs(substitutions)
    result_output_wave = result_output_wave.subs(substitutions)
    return (Quantity(result_input_wave), Quantity(result_output_wave))
