"""
Electric field due to dipole (Vector)
=====================================

"""

from sympy import Eq, Rational, pi
from symplyphysics import validate_input, validate_output, symbols, quantities

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot, VectorNorm
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

electric_field = clone_as_vector_symbol(symbols.electric_field_strength)
"""
Vector of the electric field. See :symbols:`electric_field_strength`.
"""

electric_dipole_moment = clone_as_vector_symbol(symbols.electric_dipole_moment)
"""
Vector of the :symbols:`electric_dipole_moment`.
"""

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Position vector of the point in space.  See :symbols:`distance_to_origin`.
"""

electric_field_law = Eq(
    electric_field,
    (1 / (4 * pi * quantities.vacuum_permittivity)) *
    ((3 * VectorDot(electric_dipole_moment, position_vector) / VectorNorm(position_vector)**5) *
    position_vector - electric_dipole_moment / VectorNorm(position_vector)**3),
)
"""
:laws:symbol::

:laws:latex::
"""

electric_dipole_moment_law = Eq(
    electric_dipole_moment,
    (4 * pi * quantities.vacuum_permittivity) *
    (Rational(3, 2) * VectorDot(electric_field, position_vector) * VectorNorm(position_vector) *
    position_vector - VectorNorm(position_vector)**3 * electric_field),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    electric_dipole_moment_=electric_dipole_moment,
    position_vector_=position_vector,
)
@validate_output(electric_field)
def calculate_electric_field(
    electric_dipole_moment_: QuantityCoordinateVector,
    position_vector_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = electric_field_law.rhs.subs({
        electric_dipole_moment: electric_dipole_moment_,
        position_vector: position_vector_,
    })

    return QuantityCoordinateVector.from_expr(result)


@validate_input(
    electric_field_=electric_field,
    position_vector_=position_vector,
)
@validate_output(electric_dipole_moment)
def calculate_electric_dipole_moment(
    electric_field_: QuantityCoordinateVector,
    position_vector_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = electric_dipole_moment_law.rhs.subs({
        electric_field: electric_field_,
        position_vector: position_vector_,
    })

    return QuantityCoordinateVector.from_expr(result)
