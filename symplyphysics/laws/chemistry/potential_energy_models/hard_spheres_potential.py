from sympy import Eq, Piecewise, S
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Hard spheres are widely used as model particles in the statistical mechanical theory, defined as
## impenetrable spheres that cannot overlap in space, which mimics extremely strong repulsion that
## atoms and molecules experience at very close distances

# Law: U(r) = infinity if r <= sigma else 0
## U - hard spheres potential
## r - distance between particles
## sigma - diameter of sphere

potential = Symbol("potential", units.energy)
distance = Symbol("distance", units.length)
sphere_diameter = Symbol("sphere_diameter", units.length)

law = Eq(
    potential,
    Piecewise((S.Infinity, distance <= sphere_diameter), (0, distance > sphere_diameter)),
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    distance_=distance,
    sphere_diameter_=sphere_diameter,
)
@validate_output(potential)
def calculate_potential(
    distance_: Quantity,
    sphere_diameter_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        distance: distance_,
        sphere_diameter: sphere_diameter_,
    })
    return Quantity(result)
