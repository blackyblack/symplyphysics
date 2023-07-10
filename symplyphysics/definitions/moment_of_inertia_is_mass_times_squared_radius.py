from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## If the particle is about to spin around axle, it has moment of inertia.
## Moment of inertia describes object's ability to be inertial to spinning, as object has mass to be inertial to linear movement.
# Definition: I = m * R**2
# Where:
## I is moment of inertia,
## m is mass of particle,
## R is distance to spin axle.
# Conditions
## Particle is zero-sized, rigid and uniform.

moment_of_inertia = Symbol("moment_of_inertia", units.mass * units.area)
particle_mass = Symbol("particle_mass", units.mass)
spinning_radius = Symbol("spinning_radius", units.length)

definition = Eq(moment_of_inertia, particle_mass * spinning_radius**2)

definition_units_SI = units.kilogram * units.meter**2


def print_law() -> str:
    return print_expression(definition)


@validate_input(mass_=particle_mass, radius_=spinning_radius)
@validate_output(moment_of_inertia)
def calculate_moment_of_inertia(mass_: Quantity, radius_: Quantity) -> Quantity:
    result_inertia_expr = solve(definition, moment_of_inertia, dict=True)[0][moment_of_inertia]
    result_expr = result_inertia_expr.subs({particle_mass: mass_, spinning_radius: radius_})
    return Quantity(result_expr)
