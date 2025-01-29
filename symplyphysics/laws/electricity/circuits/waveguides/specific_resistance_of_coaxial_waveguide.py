"""
Specific resistance of coaxial waveguide
========================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The specific resistance of a
coaxial waveguide depends on the radius of the outer conductor and the radius of the
inner conductor, as well as on the relative permeability of the insulator material,
frequency of signal and specific conductivity of conductor.

**Notation:**

#. :quantity_notation:`vacuum_permeability`.

..
    TODO: find link
    TODO: replace `mu_0 * mu_r` with `mu`
"""

from sympy import Eq, solve, pi, sqrt
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
    quantities,
    clone_as_symbol,
)

specific_resistance = SymbolNew("R", units.impedance / units.length)
"""
:symbols:`electrical_resistance` of coaxial waveguide per unit :symbols:`length`.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the insulator.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the signal.
"""

specific_conductance = SymbolNew("G", units.conductance / units.length)
"""
:symbols:`electrical_conductance` per unit :symbols:`length`.
"""

outer_radius = clone_as_symbol(symbols.radius, subscript="\\text{o}")
"""
:symbols:`radius` of the outer conductor.
"""

inner_radius = clone_as_symbol(symbols.radius, subscript="\\text{i}")
"""
:symbols:`radius` of the inner conductor.
"""

law = Eq(specific_resistance,
    (1 / (2 * pi)) * sqrt(angular_frequency * quantities.vacuum_permeability * relative_permeability /
    (2 * specific_conductance)) * (1 / inner_radius - 1 / outer_radius))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permeability_=relative_permeability,
    angular_frequency_=angular_frequency,
    specific_conductivity_=specific_conductance,
    outer_radius_=outer_radius,
    inner_radius_=inner_radius)
@validate_output(specific_resistance)
def calculate_specific_resistance(relative_permeability_: float, angular_frequency_: Quantity,
    specific_conductivity_: Quantity, outer_radius_: Quantity, inner_radius_: Quantity) -> Quantity:
    if outer_radius_.scale_factor <= inner_radius_.scale_factor:
        raise ValueError("The outer radius must be greater than the inner radius")
    result_expr = solve(law, specific_resistance, dict=True)[0][specific_resistance]
    result_expr = result_expr.subs({
        relative_permeability: relative_permeability_,
        angular_frequency: angular_frequency_,
        specific_conductance: specific_conductivity_,
        outer_radius: outer_radius_,
        inner_radius: inner_radius_
    })
    return Quantity(result_expr)
