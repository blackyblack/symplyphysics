"""
Specific resistance of coaxial waveguide
========================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The specific resistance of a
coaxial waveguide depends on the radius of the outer conductor and the radius of the
inner conductor, as well as on the permeability of the insulator material, frequency of
signal and specific conductivity of conductor.

..
    TODO: find link
"""

from sympy import Eq, solve, pi, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

specific_resistance = Symbol("R", units.impedance / units.length)
"""
:symbols:`electrical_resistance` of coaxial waveguide per unit :symbols:`length`.
"""

absolute_permeability = symbols.absolute_permeability
"""
:symbols:`absolute_permeability` of the insulator.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the signal.
"""

specific_conductance = Symbol("G", units.conductance / units.length)
"""
:symbols:`electrical_conductance` per unit :symbols:`length`.
"""

outer_radius = clone_as_symbol(symbols.radius, display_symbol="r_o", display_latex="r_\\text{o}")
"""
:symbols:`radius` of the outer conductor.
"""

inner_radius = clone_as_symbol(symbols.radius, display_symbol="r_i", display_latex="r_\\text{i}")
"""
:symbols:`radius` of the inner conductor.
"""

law = Eq(specific_resistance,
    (1 / (2 * pi)) * sqrt(angular_frequency * absolute_permeability /
    (2 * specific_conductance)) * (1 / inner_radius - 1 / outer_radius))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absolute_permeability_=absolute_permeability,
    angular_frequency_=angular_frequency,
    specific_conductivity_=specific_conductance,
    outer_radius_=outer_radius,
    inner_radius_=inner_radius)
@validate_output(specific_resistance)
def calculate_specific_resistance(absolute_permeability_: Quantity, angular_frequency_: Quantity,
    specific_conductivity_: Quantity, outer_radius_: Quantity, inner_radius_: Quantity) -> Quantity:
    if outer_radius_.scale_factor <= inner_radius_.scale_factor:
        raise ValueError("The outer radius must be greater than the inner radius")
    result_expr = solve(law, specific_resistance, dict=True)[0][specific_resistance]
    result_expr = result_expr.subs({
        absolute_permeability: absolute_permeability_,
        angular_frequency: angular_frequency_,
        specific_conductance: specific_conductivity_,
        outer_radius: outer_radius_,
        inner_radius: inner_radius_
    })
    return Quantity(result_expr)
