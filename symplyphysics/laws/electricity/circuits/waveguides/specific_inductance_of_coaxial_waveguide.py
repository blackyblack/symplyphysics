"""
Specific inductance of coaxial waveguide
========================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The specific inductance of a
coaxial waveguide depends on the radius of the outer conductor and the radius of the
inner conductor, as well as on the permeability of the insulator material.

..
    TODO: find link
"""

from sympy import Eq, solve, pi, ln
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

specific_inductance = Symbol("L", units.inductance / units.length)
"""
:symbols:`inductance` of the waveguide per unit :symbols:`length`.
"""

absolute_permeability = symbols.absolute_permeability
"""
:symbols:`absolute_permeability` of the insulator.
"""

outer_radius = clone_as_symbol(symbols.radius, display_symbol="r_o", display_latex="r_\\text{o}")
"""
:symbols:`radius` of the outer conductor.
"""

inner_radius = clone_as_symbol(symbols.radius, display_symbol="r_i", display_latex="r_\\text{i}")
"""
:symbols:`radius` of the inner conductor.
"""

law = Eq(specific_inductance,
    (absolute_permeability / (2 * pi)) * ln(outer_radius / inner_radius))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absolute_permeability_=absolute_permeability,
    outer_radius_=outer_radius,
    inner_radius_=inner_radius)
@validate_output(specific_inductance)
def calculate_specific_inductance(absolute_permeability_: Quantity, outer_radius_: Quantity,
    inner_radius_: Quantity) -> Quantity:
    if outer_radius_.scale_factor <= inner_radius_.scale_factor:
        raise ValueError("The outer radius must be greater than the inner radius")
    result_velocity_expr = solve(law, specific_inductance, dict=True)[0][specific_inductance]
    result_expr = result_velocity_expr.subs({
        absolute_permeability: absolute_permeability_,
        outer_radius: outer_radius_,
        inner_radius: inner_radius_
    })
    return Quantity(result_expr)
