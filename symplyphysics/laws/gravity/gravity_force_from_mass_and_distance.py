"""
Gravitational force from mass and distance
==========================================

In classical mechanics, the **gravitational force** is a fundamental attractive force that
exists between any two massive bodies. Its magnitude is proportional to the mass of the
bodies and inversely proportional to the distance squared between the bodies.

**Notation:**

#. :quantity_notation:`gravitational_constant`.

**Links:**

#. `Physics LibreTexts. Newton's Law of Universal Gravitation (13.2.1) <https://phys.libretexts.org/Workbench/PH_245_Textbook_V2/13%3A_Gravitation/13.02%3A_Newton's_Law_of_Universal_Gravitation>`__.
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    vector_magnitude,
    clone_as_symbol,
    symbols,
    quantities,
)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.points.cartesian_point import CartesianPoint
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.fields.scalar_field import ScalarField
from symplyphysics.laws.dynamics.fields import (
    conservative_force_is_gradient_of_potential_energy as gradient_law,)
from symplyphysics.laws.gravity import gravitational_potential_energy

gravitational_force = symbols.force
"""
Gravitational :symbols:`force`.
"""

first_mass = clone_as_symbol(symbols.mass, subscript="1")
"""
:symbols:`mass` of the first body.
"""

second_mass = clone_as_symbol(symbols.mass, subscript="2")
"""
:symbols:`mass` of the second body.
"""

distance_between_mass_centers = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between the centers of the bodies..
"""

law = Eq(gravitational_force,
    quantities.gravitational_constant * first_mass * second_mass / distance_between_mass_centers**2)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from the gravitational potential energy
# Condition: space must be 3-dimensional and flat

_potential = gravitational_potential_energy.law.rhs.subs({
    gravitational_potential_energy.first_mass: first_mass,
    gravitational_potential_energy.second_mass: second_mass,
    gravitational_potential_energy.distance_between_mass_centers: distance_between_mass_centers,
})


def potential_field_function(point: CartesianPoint) -> ScalarValue:
    return _potential.subs(
        distance_between_mass_centers,
        sqrt(point.x**2 + point.y**2 + point.z**2),
    )


_potential_field = ScalarField(potential_field_function)

_gravitational_force_vector = gradient_law.law(_potential_field)
_gravitational_force_derived = vector_magnitude(_gravitational_force_vector).simplify()

_x, _y, _z = _gravitational_force_vector.coordinate_system.coord_system.base_scalars()
_gravitational_force_from_law = law.rhs.subs(distance_between_mass_centers, sqrt(_x**2 + _y**2 + _z**2))

# sympy avoids oversimplifications in case of square roots without certain assumptions,
# therefore we resort to squaring both sides to make it work
assert expr_equals(_gravitational_force_derived**2, _gravitational_force_from_law**2)


@validate_input(first_object_mass_=first_mass,
    second_object_mass_=second_mass,
    distance_between_objects_=distance_between_mass_centers)
@validate_output(gravitational_force)
def calculate_force(first_object_mass_: Quantity, second_object_mass_: Quantity,
    distance_between_objects_: Quantity) -> Quantity:
    result_force_expr = solve(law, gravitational_force, dict=True)[0][gravitational_force]
    result_expr = result_force_expr.subs({
        first_mass: first_object_mass_,
        second_mass: second_object_mass_,
        distance_between_mass_centers: distance_between_objects_
    })
    return Quantity(result_expr)
