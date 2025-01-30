"""
Attenuation coefficient in metal
================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The specific resistance of a
coaxial waveguide depends on the diameter of the outer conductor and the diameter of the
inner conductor, as well as on the relative permeability and the relative permittivity
of the insulator material, the surface resistance of the outer conductor and the surface
resistance of the inner conductor.

..
    TODO: find link
"""

from sympy import Eq, solve, pi, sqrt, ln
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

attenuation_coefficient = symbols.attenuation_coefficient
"""
:symbols:`attenuation_coefficient` in metal.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the insulator.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the insulator.
"""

outer_surface_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="\\text{o}")
"""
Surface :symbols:`electrical_resistance` of the outer conductor.
"""

inner_surface_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="\\text{i}")
"""
Surface :symbols:`electrical_resistance` of the inner conductor.
"""

outer_diameter = clone_as_symbol(symbols.diameter, subscript="\\text{o}")
"""
:symbols:`diameter` of the outer conductor.
"""

inner_diameter = clone_as_symbol(symbols.diameter, subscript="\\text{i}")
"""
:symbols:`diameter` of the inner conductor.
"""

resistance = Quantity(420 * units.ohm, display_symbol="R_0")
"""
Constant equal to :math:`420 \\Omega`.
"""

law = Eq(
    attenuation_coefficient,
    sqrt(relative_permittivity / relative_permeability) *
    ((inner_surface_resistance / inner_diameter) + (outer_surface_resistance / outer_diameter)) /
    (pi * resistance * ln(outer_diameter / inner_diameter)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability,
    surface_resistance_outer_=outer_surface_resistance,
    inner_surface_resistance_=inner_surface_resistance,
    outer_diameter_=outer_diameter,
    inner_diameter_=inner_diameter)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(relative_permittivity_: float, relative_permeability_: float,
    surface_resistance_outer_: Quantity, inner_surface_resistance_: Quantity,
    outer_diameter_: Quantity, inner_diameter_: Quantity) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    if outer_diameter_.scale_factor <= inner_diameter_.scale_factor:
        raise ValueError("The outer diameter must be greater than the inner diameter")
    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
        outer_surface_resistance: surface_resistance_outer_,
        inner_surface_resistance: inner_surface_resistance_,
        outer_diameter: outer_diameter_,
        inner_diameter: inner_diameter_
    })
    return Quantity(result_expr)
