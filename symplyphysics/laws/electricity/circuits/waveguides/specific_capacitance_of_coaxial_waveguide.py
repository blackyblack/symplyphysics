"""
Specific capacitance of coaxial waveguide
=========================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The specific capacitance of a
coaxial waveguide depends on the radius of the outer conductor and the radius of the
inner conductor, as well as on the permittivity of the insulator material.

..
    TODO: find link
"""

from sympy import Eq, solve, pi, ln
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

specific_capacitance = SymbolNew("C", units.capacitance / units.length)
"""
:symbols:`capacitance` of the waveguide per unit :symbols:`length`.
"""

absolute_permittivity = symbols.absolute_permittivity
"""
:symbols:`absolute_permittivity` of the insulator.
"""

outer_radius = clone_as_symbol(symbols.radius, display_symbol="r_o", display_latex="r_\\text{o}")
"""
:symbols:`radius` of the outer conductor.
"""

inner_radius = clone_as_symbol(symbols.radius, display_symbol="r_i", display_latex="r_\\text{i}")
"""
:symbols:`radius` of the inner conductor.
"""

law = Eq(specific_capacitance,
    (2 * pi * absolute_permittivity) / ln(outer_radius / inner_radius))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absolute_permittivity_=absolute_permittivity,
    outer_radius_=outer_radius,
    inner_radius_=inner_radius)
@validate_output(specific_capacitance)
def calculate_specific_capacitance(absolute_permittivity_: Quantity, outer_radius_: Quantity,
    inner_radius_: Quantity) -> Quantity:
    if outer_radius_.scale_factor <= inner_radius_.scale_factor:
        raise ValueError("The outer radius must be greater than the inner radius")
    result_velocity_expr = solve(law, specific_capacitance, dict=True)[0][specific_capacitance]
    result_expr = result_velocity_expr.subs({
        absolute_permittivity: absolute_permittivity_,
        outer_radius: outer_radius_,
        inner_radius: inner_radius_
    })
    return Quantity(result_expr)
