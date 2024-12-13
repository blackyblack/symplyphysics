from sympy import Eq, solve, pi
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The effective cross section is a physical quantity characterizing the probability of transition of a system of
## two interacting particles to a certain final state, a quantitative characteristic of the acts of collision of
## particles of a stream hitting a target with target particles. The effective cross-section has the dimension of the area.

## Law is: g =  pi * d^2, where
## g - the cross-sectional area of the interaction of particles,
## d - the distance of the greatest convergence of two colliding particles

# Notes
## 1. `d = 2 * r` when the particles are spheres of radius `r`

# Links: Wikipedia, second formula <https://en.wikipedia.org/wiki/Cross_section_(physics)#Collision_among_gas_particles>

# Links: Wikipedia, second formula <https://en.wikipedia.org/wiki/Cross_section_(physics)#Collision_among_gas_particles>

# FIXME: shouldn't it be `(2*r)^2` as per the formula in the link?

cross_sectional_area_of_interaction = Symbol("cross_sectional_area_of_interaction", units.area)

distance_of_convergence_of_two_particles = Symbol("distance_of_convergence_of_two_particles",
    units.length)

law = Eq(cross_sectional_area_of_interaction, pi * distance_of_convergence_of_two_particles**2)


@validate_input(distance_of_convergence_of_two_particles_=distance_of_convergence_of_two_particles)
@validate_output(cross_sectional_area_of_interaction)
def calculate_cross_sectional_area_of_interaction(
        distance_of_convergence_of_two_particles_: Quantity) -> Quantity:
    result_expr = solve(law, cross_sectional_area_of_interaction,
        dict=True)[0][cross_sectional_area_of_interaction]
    result_expr = result_expr.subs({
        distance_of_convergence_of_two_particles: distance_of_convergence_of_two_particles_,
    })
    return Quantity(result_expr)
