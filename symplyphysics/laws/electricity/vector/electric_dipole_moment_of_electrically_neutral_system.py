"""
Electric dipole moment of electrically neutral system
=====================================================

The electric dipole moment of an electrically neutral system of (point) charges can be found as
the sum of the products of the values and the position vectors of the charges that compose the
system.

**Notes:**

#. This law can be seen as the definition of electric dipole moment.

#. The value of the electric dipole moment for such a system is independent of the choice of the
   origin of the coordinate frame (i.e. it is translationally invariant).

..
    TODO: derive this property of the electric dipole moment vector

**Conditions:**

#. The system is electrically neutral.

**Links:**

#. `Wikipedia, derivable from the third equation <https://en.wikipedia.org/wiki/Electric_dipole_moment#Expression_(general_case)>`__.
"""

from typing import Sequence, Optional
from sympy import Eq, Idx
from symplyphysics import (Quantity, validate_input, validate_output, symbols, global_index,
    IndexedSum, assert_equal)
from symplyphysics.core.symbols.symbols import clone_as_indexed

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, clone_as_indexed_vector

electric_dipole_moment = clone_as_vector_symbol(symbols.electric_dipole_moment)
"""
Vector of the :symbols:`electric_dipole_moment` of the system of charges.
"""

charge = clone_as_indexed(symbols.charge)
"""
Value of the :math:`i`-th point charge.
"""

position_vector = clone_as_indexed_vector(symbols.distance_to_origin)
"""
Position vector of the :math:`i`-th point charge. See :symbols:`distance_to_origin`.
"""

law = Eq(
    electric_dipole_moment,
    IndexedSum(charge[global_index] * position_vector[global_index], global_index),
)
"""
:laws:symbol::

:laws:latex::
"""

electric_neutrality_condition = Eq(
    IndexedSum(charge[global_index], global_index),
    0,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    charges_=charge,
    position_vectors_=position_vector,
)
@validate_output(electric_dipole_moment)
def calculate_electric_dipole_moment(
    charges_: Sequence[Quantity],
    position_vectors_: Sequence[QuantityCoordinateVector],
    absolute_tolerance_: Optional[float] = None,
) -> QuantityCoordinateVector:
    n = len(charges_)

    if n != len(position_vectors_):
        raise ValueError("Number of charges must match the number of position vectors.")

    # See `electric_neutrality_condition`
    assert_equal(sum(charges_), 0, absolute_tolerance=absolute_tolerance_)

    local_index = Idx("i", (1, n))

    result = law.rhs.subs(global_index, local_index).doit()

    for i, (charge_, position_vector_) in enumerate(zip(charges_, position_vectors_), start=1):
        result = result.subs({
            charge[i]: charge_,
            position_vector[i]: position_vector_,
        })

    return QuantityCoordinateVector.from_expr(result)


# UNIQUE_LAW_ID: 536
