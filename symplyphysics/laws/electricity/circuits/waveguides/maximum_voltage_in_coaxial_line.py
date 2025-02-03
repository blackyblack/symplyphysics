"""
Maximum voltage in coaxial line
===============================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. An electrical breakdown is a
phenomenon of a sharp increase in current that occurs when the field intensity is higher
than the critical one — dielectric breakdown intensity.

..
    TODO: find link
"""

from sympy import Eq, solve, ln
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

maximum_voltage = symbols.voltage
"""
Maximum :symbols:`voltage` between the central conductor and the outer conductor.
"""

breakdown_electric_field = symbols.electric_field_strength
"""
:symbols:`electric_field_strength` of dielectric breakdown.
"""

outer_diameter = clone_as_symbol(symbols.diameter, display_symbol="d_o", display_latex="d_\\text{o}")
"""
:symbols:`diameter` of the outer conductor.
"""

inner_diameter = clone_as_symbol(symbols.diameter, display_symbol="d_i", display_latex="d_\\text{i}")
"""
:symbols:`diameter` of the inner conductor.
"""

law = Eq(maximum_voltage,
    breakdown_electric_field * outer_diameter * ln(outer_diameter / inner_diameter) / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(breakdown_intensity_=breakdown_electric_field,
    outer_diameter_=outer_diameter,
    inner_diameter_=inner_diameter)
@validate_output(maximum_voltage)
def calculate_maximum_voltage(breakdown_intensity_: Quantity, outer_diameter_: Quantity,
    inner_diameter_: Quantity) -> Quantity:
    if outer_diameter_.scale_factor <= inner_diameter_.scale_factor:
        raise ValueError("The outer diameter must be greater than the inner diameter")
    result_velocity_expr = solve(law, maximum_voltage, dict=True)[0][maximum_voltage]
    result_expr = result_velocity_expr.subs({
        breakdown_electric_field: breakdown_intensity_,
        outer_diameter: outer_diameter_,
        inner_diameter: inner_diameter_
    })
    return Quantity(result_expr)
