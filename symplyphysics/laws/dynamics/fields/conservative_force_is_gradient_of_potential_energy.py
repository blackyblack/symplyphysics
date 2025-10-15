"""
Conservative force is gradient of potential energy
==================================================

A conservative force is a such a force, the total work of which in moving a particle between two
points is independent of the path taken. Alternative definition states that if a particle travels
in a closed loop, the total work done by a conservative force is zero.

**Conditions:**

#. Force is conservative. Mathematically, this can be expressed as
   :math:`\\text{curl} \\, {\\vec F} \\! \\left( \\vec r \\right) \\equiv 0`, i.e. the force field must
   be irrotational.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Conservative_force>`__.
"""

from sympy import Eq

from symplyphysics import symbols, clone_as_function

from symplyphysics.core.experimental.vectors import (clone_as_vector_symbol,
    clone_as_vector_function)
from symplyphysics.core.experimental.operators import VectorGradient

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Position vector of a point in space. See :symbols:`distance_to_origin`.
"""

force = clone_as_vector_function(symbols.force, (position_vector,))
"""
Vector field of the conservative force as a function of the :attr`position_vector`. See
:symbols:`force`.
"""

potential_energy = clone_as_function(symbols.potential_energy, (position_vector,))
"""
Scalar field of the force's potential as a function of the :attr`position_vector`. See
:symbols:`potential_energy`.
"""

law = Eq(
    force(position_vector),
    -1 * VectorGradient(potential_energy(position_vector), evaluate=False),
)
"""
:laws:symbol::

:laws:latex::
"""


# UNIQUE_LAW_ID: 219
