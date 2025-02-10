"""
Resonant frequency of rectangular resonator
===========================================

A rectangular resonator consists of metal walls and a material filling it.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Links:**

#. `Mahatma Gandhi Central University, formula 17 on page 12 (PDF file) <https://mgcub.ac.in/pdf/material/20200427120209be39415ad1.pdf>`__.
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    dimensionless,
    symbols,
    quantities,
    clone_as_symbol,
)

resonant_frequency = clone_as_symbol(symbols.temporal_frequency, display_symbol="f_r", display_latex="f_\\text{r}")
"""
Resonant :symbols:`temporal_frequency` of the resonator.
"""

first_index = clone_as_symbol(symbols.positive_number, display_symbol="m", display_latex="m")
"""
Index that changes along the :attr:`~width` of the resonator, see
:symbols:`positive_number`.
"""

second_index = clone_as_symbol(symbols.positive_number, display_symbol="n", display_latex="n")
"""
Index that changes along the :attr:`~height` of the resonator, see
:symbols:`positive_number`.
"""

third_index = clone_as_symbol(symbols.positive_number, display_symbol="p", display_latex="p")
"""
Index that changes along the :attr:`~length` of the resonator, see
:symbols:`positive_number`.
"""

width = clone_as_symbol(symbols.length, subscript="1")
"""
:symbols:`length` of the resonator along the axis perpendicular to the axis of wave
propagation and to :attr:`~height`.
"""

height = clone_as_symbol(symbols.length, subscript="2")
"""
:symbols:`length` of the resonator along the axis perpendicular to the axis of wave
propagation and to :attr:`~width`.
"""

length = clone_as_symbol(symbols.length, subscript="3")
"""
:symbols:`length` of the resonator along the axis of wave propagation.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the medium filling the resonator.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the medium filling the resonator.
"""

law = Eq(
    resonant_frequency,
    (quantities.speed_of_light / (2 * sqrt(relative_permittivity * relative_permeability)))
    * sqrt((first_index / width)**2 + (second_index / height)**2 + (third_index / length)**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(indexes_=dimensionless,
    resonator_dimensions_=units.length,
    relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability)
@validate_output(resonant_frequency)
def calculate_resonant_frequency(indexes_: tuple[float, float,
    float], resonator_dimensions_: tuple[Quantity, Quantity, Quantity],
    relative_permittivity_: Quantity, relative_permeability_: Quantity) -> Quantity:
    first_index_, second_index_, third_index_ = indexes_
    resonator_width_, resonator_height_, resonator_length_ = resonator_dimensions_
    result_expr = solve(law, resonant_frequency, dict=True)[0][resonant_frequency]
    result_expr = result_expr.subs({
        first_index: first_index_,
        second_index: second_index_,
        third_index: third_index_,
        width: resonator_width_,
        height: resonator_height_,
        length: resonator_length_,
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
    })
    return Quantity(result_expr)
