"""
Quality factor of empty rectangular resonator for traverse electric waves
=========================================================================

A rectangular resonator consists of metal walls and a material filling it. In the case
when the resonator is empty, its quality factor depends only on the losses in the metal
walls of the resonator.

**Notation:**

#. :quantity_notation:`vacuum_permeability`.

**Conditions:**

#. Applies to transverse electric waves with indices :math:`m = p = 1, n = 0`.

..
    TODO: find link
    TODO: replace `mu_0 * mu_r` with `mu`
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
    quantities,
)

quality_factor = symbols.quality_factor
"""
:symbols:`quality_factor` of the resonator.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the current.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the resonator walls.
"""

surface_resistance = clone_as_symbol(symbols.electrical_resistance, display_symbol="R_s", display_latex="R_\\text{s}")
"""
:symbols:`electrical_resistance` of the resonator walls.
"""

length = clone_as_symbol(symbols.length, subscript="1")
"""
:symbols:`length` of the resonator along the axis of wave propagation.
"""

width = clone_as_symbol(symbols.length, subscript="2")
"""
:symbols:`length` of the resonator along the axis perpendicular to the axis of wave
propagation and to :attr:`~height`.
"""

height = clone_as_symbol(symbols.length, subscript="3")
"""
:symbols:`length` of the resonator along the axis perpendicular to the axis of wave
propagation and to :attr:`~width`.
"""

law = Eq(
    quality_factor,
    angular_frequency * quantities.vacuum_permeability * relative_permeability * height *
    width * length * (width**2 + length**2) /
    (2 * surface_resistance * (width**3 * (length + 2 * height) + length**3 * (width + 2 * height))))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(angular_frequency_=angular_frequency,
    relative_permeability_=relative_permeability,
    surface_resistance_=surface_resistance,
    resonator_dimensions_=units.length)
@validate_output(quality_factor)
def calculate_quality_factor(angular_frequency_: Quantity, relative_permeability_: float,
    surface_resistance_: Quantity, resonator_dimensions_: tuple[Quantity, Quantity,
    Quantity]) -> float:
    resonator_width_, resonator_height_, resonator_length_ = resonator_dimensions_
    result_expr = solve(law, quality_factor, dict=True)[0][quality_factor]
    result_expr = result_expr.subs({
        angular_frequency: angular_frequency_,
        relative_permeability: relative_permeability_,
        surface_resistance: surface_resistance_,
        width: resonator_width_,
        height: resonator_height_,
        length: resonator_length_
    })
    return convert_to_float(result_expr)
