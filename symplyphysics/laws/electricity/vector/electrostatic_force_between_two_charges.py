"""
Electrostatic force between two charges
=======================================

Also known as **Coulomb's law**, it is an experimental law that calculates the amount
of force between two electrically charged particles at rest.

Also see the :ref:`scalar law <Electrostatic force via charges and distance>`.

**Notes:**

#. If a given charge is in the vicinity of a system of point charges, then the net law can be
   found via the :ref:`principle of superposition <Superposition of forces is sum (Vector)>`.

**Notation:**

#. :quantity_notation:`vacuum_permittivity`.

**Conditions:**

#. The charges are small.

#. The charges are at rest.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Coulomb%27s_law#Mathematical_form>`__.
"""

from sympy import Eq, pi, sign, sqrt, Rational
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol,
    quantities)

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorNorm
from symplyphysics.core.experimental.solvers import solve_for_vector

force = clone_as_vector_symbol(symbols.force, subscript="12")
"""
Vector of the electrostatic :symbols:`force` experienced by the :attr:`~first_charge` in the
vicinity of the :attr:`~second_charge` in vacuum.
"""

first_charge = clone_as_symbol(symbols.charge, subscript="1")
"""
Value of the first point :symbols:`charge`.
"""

second_charge = clone_as_symbol(symbols.charge, subscript="2")
"""
Value of the second point :symbols:`charge`.
"""

position_vector = clone_as_vector_symbol(symbols.euclidean_distance, subscript="21")
"""
Position vector drawn from the :attr:`~second_charge` to the :attr:`~first_charge`.
"""

force_law = Eq(
    force,
    ((first_charge * second_charge) / (4 * pi * quantities.vacuum_permittivity)) *
    (position_vector / VectorNorm(position_vector)**3),
)
"""
:laws:symbol::

:laws:latex::
"""

position_vector_law = Eq(
    position_vector,
    sign(first_charge * second_charge) *
    sqrt(abs(first_charge * second_charge) /
    (4 * pi * quantities.vacuum_permittivity)) * (force / VectorNorm(force)**Rational(3, 2)),
)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive `position_vector_law` from `force_law`


@validate_input(
    first_charge_=first_charge,
    second_charge_=second_charge,
    position_vector_=position_vector,
)
@validate_output(force)
def calculate_force(
    first_charge_: Quantity,
    second_charge_: Quantity,
    position_vector_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = solve_for_vector(force_law, force).subs({
        first_charge: first_charge_,
        second_charge: second_charge_,
        position_vector: position_vector_,
    })

    return QuantityCoordinateVector.from_expr(result)


@validate_input(
    first_charge_=first_charge,
    second_charge_=second_charge,
    force_=force,
)
@validate_output(position_vector)
def calculate_position_vector(
    first_charge_: Quantity,
    second_charge_: Quantity,
    force_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = position_vector_law.rhs.subs({
        first_charge: first_charge_,
        second_charge: second_charge_,
        force: force_,
    })

    return QuantityCoordinateVector.from_expr(result)
