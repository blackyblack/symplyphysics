from sympy import Eq, solve, sqrt
from sympy.physics.units import gravitational_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    vector_magnitude,
    clone_as_symbol,
    symbols,
)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.points.cartesian_point import CartesianPoint
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.fields.scalar_field import ScalarField
from symplyphysics.laws.dynamics.fields import (
    conservative_force_is_gradient_of_potential_energy as gradient_law,)
from symplyphysics.laws.gravity import gravitational_potential_energy

# Description
## Every object generates gravity field around it. Any other object in this field is pulled toward generator.
## Gravitational force between two object if proportional to masses of objects and counter-proportional to distance between their mass centers.
## Law: F = G * m1 * m2 / R**2
## Where:
## F - gravitational force
## m1 and m2 - masses of objects
## R - distance between mass centers of objects
## G - gravitational constant

gravitational_force = symbols.force
first_mass = clone_as_symbol(symbols.mass, display_symbol="m_1", display_latex="m_1")
second_mass = clone_as_symbol(symbols.mass, display_symbol="m_2", display_latex="m_2")
distance_between_mass_centers = Symbol("distance_between_mass_centers", units.length)

law = Eq(gravitational_force,
    gravitational_constant * first_mass * second_mass / distance_between_mass_centers**2)

# Derive law from the gravitational potential energy
# Condition: space must be 3-dimensional and flat

potential = gravitational_potential_energy.law.rhs.subs({
    gravitational_potential_energy.first_mass: first_mass,
    gravitational_potential_energy.second_mass: second_mass,
    gravitational_potential_energy.distance_between_mass_centers: distance_between_mass_centers,
})


def potential_field_function(point: CartesianPoint) -> ScalarValue:
    return potential.subs(
        distance_between_mass_centers,
        sqrt(point.x**2 + point.y**2 + point.z**2),
    )


potential_field = ScalarField(potential_field_function)

gravitational_force_vector = gradient_law.law(potential_field)
gravitational_force_derived = vector_magnitude(gravitational_force_vector).simplify()

x, y, z = gravitational_force_vector.coordinate_system.coord_system.base_scalars()
gravitational_force_from_law = law.rhs.subs(distance_between_mass_centers, sqrt(x**2 + y**2 + z**2))

# sympy avoids oversimplifications in case of square roots without certain assumptions,
# therefore we resort to squaring both sides to make it work
assert expr_equals(gravitational_force_derived**2, gravitational_force_from_law**2)


def print_law() -> str:
    return print_expression(law)


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
