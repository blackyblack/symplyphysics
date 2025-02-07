"""
Resonant frequency of rectangular resonator
===========================================

A rectangular resonator consists of metal walls and a material filling it. The indices
show the number of half-waves that fit along the length, width, height of the resonator,
respectively.

**Notation:**

#. :quantity_notation:`speed_of_light`.

..
    TODO: find link
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

resonant_frequency = symbols.temporal_frequency
"""
Resonant :symbols:`temporal_frequency` of the resonator.
"""

first_index = clone_as_symbol(symbols.positive_number, subscript="1")
"""
First index, see :symbols:`positive_number`.
"""

second_index = clone_as_symbol(symbols.positive_number, subscript="2")
"""
Second index, see :symbols:`positive_number`.
"""

third_index = clone_as_symbol(symbols.positive_number, subscript="3")
"""
Third index, see :symbols:`positive_number`.
"""

resonator_length = clone_as_symbol(symbols.length, subscript="1")
"""
:symbols:`length` of the resonator along the first dimension.
"""

resonator_width = clone_as_symbol(symbols.length, subscript="2")
"""
:symbols:`length` of the resonator along the second dimension.
"""

resonator_height = clone_as_symbol(symbols.length, subscript="3")
"""
:symbols:`length` of the resonator along the third dimension.
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
    quantities.speed_of_light 
    * sqrt((first_index / resonator_length)**2 + (second_index / resonator_width)**2 + (third_index / resonator_height)**2) 
    / (2 * sqrt(relative_permittivity * relative_permeability)))
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
        resonator_width: resonator_width_,
        resonator_height: resonator_height_,
        resonator_length: resonator_length_,
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
    })
    return Quantity(result_expr)
