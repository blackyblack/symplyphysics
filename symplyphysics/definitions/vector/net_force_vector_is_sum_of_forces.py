r"""
Net force vector is sum of forces
=================================

The net force exerted on an object is equal to the vector sum of all the forces exerted on it.

**Notation:**

#. :math:`\sum_i x_i` (:code:`Sum(x_i, i)`) denotes a sum of :math:`x_i` over the index :math:`i`.

**Links:**

#. `Wikipedia, general principle <https://en.wikipedia.org/wiki/Superposition_principle>`__.
"""

from typing import Sequence

from sympy import Eq, Idx

from symplyphysics import symbols, validate_input, validate_output, IndexedSum, global_index

from symplyphysics.core.experimental.solvers import solve_for_vector
from symplyphysics.core.experimental.vectors import (clone_as_vector_symbol,
    clone_as_indexed_vector)
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

net_force = clone_as_vector_symbol(symbols.force)
"""
Vector of the net :symbols:`force` exerted on the body.
"""

force = clone_as_indexed_vector(symbols.force)
"""
Vector of an individual :symbols:`force`.
"""

law = Eq(net_force, IndexedSum(force[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(forces_=force)
@validate_output(net_force)
def calculate_resultant_force(
        forces_: Sequence[QuantityCoordinateVector]) -> QuantityCoordinateVector:
    local_index = Idx("i", (1, len(forces_)))
    net_force_value = solve_for_vector(law, net_force).subs(global_index, local_index).doit()

    for i, force_i in enumerate(forces_, start=1):
        net_force_value = net_force_value.subs(force[i], force_i)

    return QuantityCoordinateVector.from_expr(net_force_value)
