"""
Center of mass for a system of particles
========================================

The center of mass (com) of a system of particles is a unique point at any given time where
the sum of weighted relative positions of the distributed mass is zero.

**Links:**

#. `Wikipedia, second formula <https://en.wikipedia.org/wiki/Center_of_mass#A_system_of_particles>`__.
"""

from typing import Sequence
from sympy import Eq, Idx

from symplyphysics import (validate_input, validate_output, Quantity, symbols, global_index,
    IndexedSum, units)
from symplyphysics.core.symbols.symbols import clone_as_indexed

from symplyphysics.core.experimental.vectors import clone_as_indexed_vector, VectorSymbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

center_of_mass = VectorSymbol("r_com", units.length, display_latex="{\\vec r}_\\text{COM}")
"""
Vector of the system's center of mass (COM).
"""

position_vector = clone_as_indexed_vector(symbols.distance_to_origin)
"""
Position vector of the :math:`i`-th body. See :symbols:`distance_to_origin`.
"""

mass = clone_as_indexed(symbols.mass)
"""
:symbols:`mass` of the :math:`i`-th body.
"""

law = Eq(
    center_of_mass,
    IndexedSum(mass[global_index] * position_vector[global_index], global_index) /
    IndexedSum(mass[global_index], global_index),
)
"""
:laws:symbol::

..
    For the Latex code printer:
    TODO: fix indexed vector symbols
    TODO: add parenthesis around IndexedSum when it is the base of an exponent

Latex:
    .. math::
        {\\vec r}_\\text{COM} = \\frac{\\sum_i m_i {\\vec r}_i}{\\sum_i m_i}
"""


@validate_input(masses_=mass, position_vectors_=position_vector)
@validate_output(center_of_mass)
def calculate_center_of_mass(
    masses_: Sequence[Quantity],
    position_vectors_: Sequence[QuantityCoordinateVector],
) -> QuantityCoordinateVector:
    if len(masses_) != len(position_vectors_):
        raise ValueError("Mass and position arrays should have the same lengths")

    local_index = Idx("i", (1, len(masses_)))
    result = law.rhs.subs(global_index, local_index).doit()

    for index_, (mass_, position_vector_) in enumerate(zip(masses_, position_vectors_), start=1):
        result = result.subs({
            position_vector[index_]: position_vector_,
            mass[index_]: mass_,
        })

    return QuantityCoordinateVector.from_expr(result)


# UNIQUE_LAW_ID: 443
