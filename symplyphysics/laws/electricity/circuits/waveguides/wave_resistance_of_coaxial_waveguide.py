"""
Wave resistance of coaxial waveguide
====================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The wave resistance of a coaxial
waveguide depends on the radius of the outer conductor and the radius of the inner
conductor, as well as on the relative permittivity and the relative permeability of the
insulator material.

**Notation:**

#. :quantity_notation:`vacuum_permeability`.
#. :quantity_notation:`vacuum_permittivity`.

"""

from sympy import Eq, solve, pi, sqrt, ln
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.quantities import vacuum_permeability, vacuum_permittivity

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the wave in the waveguide.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the insulator.
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

law = Eq(resistance, (1 / (2 * pi)) * sqrt(vacuum_permeability * relative_permeability /
    (vacuum_permittivity * relative_permittivity)) * ln(outer_radius / inner_radius))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability,
    outer_radius_=outer_radius,
    inner_radius_=inner_radius)
@validate_output(resistance)
def calculate_wave_resistance(relative_permittivity_: float, relative_permeability_: float,
    outer_radius_: Quantity, inner_radius_: Quantity) -> Quantity:
    if outer_radius_.scale_factor <= inner_radius_.scale_factor:
        raise ValueError("The outer radius must be greater than the inner radius")
    result_velocity_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_velocity_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
        outer_radius: outer_radius_,
        inner_radius: inner_radius_
    })
    return Quantity(result_expr)
