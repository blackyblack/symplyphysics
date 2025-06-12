"""
Conservative force is gradient of potential energy
==================================================

A conservative force is a such a force, the total work of which in moving a particle between two
points is independent of the path taken. Alternative definition states that if a particle travels
in a closed loop, the total work done by a conservative force is zero.

**Conditions:**

#. Force is conservative. Mathematically, this can be expressed as
   :math:`\\text{curl} \, {\\vec F} \\equiv 0`, i.e. the force field must be irrotational.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Conservative_force>`__.
"""

from sympy import Eq

from symplyphysics import symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.operators import VectorGradient

force = clone_as_vector_symbol(symbols.force)
"""
Vector field of the conservative force. See :symbols:`force`.
"""

potential_energy = symbols.potential_energy
"""
Scalar field of the force's potential. See :symbols:`potential_energy`.
"""

law = Eq(force, -1 * VectorGradient(potential_energy, evaluate=False))
"""
..
    NOTE: code printers have not been implemented yet for `VectorGradient`

:code:`F = -grad(U)`

Latex:
    .. math::
        \\vec F = - \\text{grad} \, U
"""
