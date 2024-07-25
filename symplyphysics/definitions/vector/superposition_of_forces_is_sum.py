r"""
Superposition of forces is sum (Vector)
=======================================

The net force exerted on an object is equal to the vector sum of all the forces exerted on it.

**Notation:**

#. :math:`\sum_i x_i` (:code:`Sum(x_i, i)`) denotes a sum of :math:`x_i` over the index :math:`i`.
"""

from typing import Iterable
from symplyphysics import (units, validate_input, validate_output)
from symplyphysics.core.vectors.arithmetics import add_cartesian_vectors
from symplyphysics.core.vectors.vectors import QuantityVector, Vector


def superposition_law(forces_: Iterable[Vector]) -> Vector:
    r"""
    The net force vector is the sum of individual force vectors.

    Law:
        :code:`F = Sum(F_i, i)`

    Latex:
        .. math::
            \vec F = \sum_i {\vec F}_i

    :param forces\_: sequence of individual :attr:`~symplyphysics.symbols.dynamics.force` vectors.

        Symbol: :code:`F_i`

        Latex: :math:`{\vec F}_i`

        Dimension: :math:`\mathsf{M} \mathsf{L} \mathsf{T}^{-2}`

    :return: net :attr:`~symplyphysics.symbols.dynamics.force` exerted on the object.

        Symbol: :code:`F`

        Latex: :math:`\vec F`

        Dimension: :math:`\mathsf{M} \mathsf{L} \mathsf{T}^{-2}`
    """

    result = Vector([0, 0, 0], next(iter(forces_)).coordinate_system)

    for force_ in forces_:
        result = add_cartesian_vectors(result, force_)

    return result


@validate_input(forces_=units.force)
@validate_output(units.force)
def calculate_resultant_force(forces_: Iterable[QuantityVector]) -> QuantityVector:
    forces_base_vectors = [f.to_base_vector() for f in forces_]
    if not forces_base_vectors:
        raise ValueError("At least one force should be present")
    result_force_vector = superposition_law(forces_base_vectors)
    return QuantityVector.from_base_vector(result_force_vector)
