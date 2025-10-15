"""
Electric field is force over test charge (Vector)
=================================================

Electric field at a point in space can be found by placing there a test charge and measuring
the electrostatic force that is applied to it.

**Links:**

#. `Wikipedia, last formula in paragraph <https://en.wikipedia.org/wiki/Electric_field#Electrostatics>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.solvers import solve_for_vector

test_charge = clone_as_symbol(symbols.charge, subscript="0")
"""
Value of the test :symbols:`charge`.
"""

electric_field = clone_as_vector_symbol(symbols.electric_field_strength)
"""
Vector of the electric field. See :symbols:`electric_field_strength`.
"""

electrostatic_force = clone_as_vector_symbol(symbols.force)
"""
Vector of the electrostatic :symbols:`force`.
"""

law = Eq(electric_field, electrostatic_force / test_charge)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(electrostatic_force_=electrostatic_force, test_charge_=test_charge)
@validate_output(electric_field)
def calculate_electric_field(
    electrostatic_force_: QuantityCoordinateVector,
    test_charge_: Quantity,
) -> QuantityCoordinateVector:
    expr = solve_for_vector(law, electric_field)
    value = expr.subs({
        electrostatic_force: electrostatic_force_,
        test_charge: test_charge_,
    })

    return QuantityCoordinateVector.from_expr(value)


@validate_input(electric_field_=electric_field, test_charge_=test_charge)
@validate_output(electrostatic_force)
def calculate_electrostatic_force(
    electric_field_: QuantityCoordinateVector,
    test_charge_: Quantity,
) -> QuantityCoordinateVector:
    expr = solve_for_vector(law, electrostatic_force)
    value = expr.subs({
        electric_field: electric_field_,
        test_charge: test_charge_,
    })

    return QuantityCoordinateVector.from_expr(value)


# UNIQUE_LAW_ID: 529
