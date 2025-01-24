"""
Specific inductance of coaxial waveguide
========================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The specific inductance of a
coaxial waveguide depends on the radius of the outer conductor and the radius of the
inner conductor, as well as on the relative permeability of the insulator material.

**Notation:**

#. :quantity_notation:`vacuum_permeability`.

..
    TODO: find link
    TODO: replace `mu_0 * mu_r` with `mu`
"""

from sympy import Eq, solve, pi, ln
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    quantities,
    symbols,
    clone_as_symbol,
)

specific_inductance = SymbolNew("L", units.inductance / units.length)
"""
:symbols:`inductance` of the waveguide per unit :symbols:`length`.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the insulator.
"""

outer_radius = clone_as_symbol(symbols.radius, subscript="\\text{o}")
"""
:symbols:`radius` of the outer conductor.
"""

inner_radius = clone_as_symbol(symbols.radius, subscript="\\text{i}")
"""
:symbols:`radius` of the inner conductor.
"""

law = Eq(specific_inductance,
    (quantities.vacuum_permeability * relative_permeability / (2 * pi)) * ln(outer_radius / inner_radius))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permeability_=relative_permeability,
    outer_radius_=outer_radius,
    inner_radius_=inner_radius)
@validate_output(specific_inductance)
def calculate_specific_inductance(relative_permeability_: float, outer_radius_: Quantity,
    inner_radius_: Quantity) -> Quantity:
    if outer_radius_.scale_factor <= inner_radius_.scale_factor:
        raise ValueError("The outer radius must be greater than the inner radius")
    result_velocity_expr = solve(law, specific_inductance, dict=True)[0][specific_inductance]
    result_expr = result_velocity_expr.subs({
        relative_permeability: relative_permeability_,
        outer_radius: outer_radius_,
        inner_radius: inner_radius_
    })
    return Quantity(result_expr)
