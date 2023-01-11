from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## If the particle is about to spin around axle, it has moment of inertia.
## Moment of inertia describes object's ability to be inertial to spinning, as object has mass to be inertial to linear movement.
## Definition: I = m * R**2
## Where I is moment of inertia,
## m is mass of particle,
## R is distance to spin axle.
# Conditions
## Particle is zero-sized, rigid and uniform.

moment_of_inertia, particle_mass, spinning_radius = symbols('moment_of_inertia particle_mass spinning_radius')
definition = Eq(moment_of_inertia, particle_mass * spinning_radius**2)

definition_dimension_SI = units.kilogram * units.meter**2

def print():
    return pretty(definition, use_unicode=False)

@validate_input(mass_=units.mass, radius_=units.length)
@validate_output(units.mass * units.length**2)
def calculate_moment_of_inertia(mass_: Quantity, radius_: Quantity) -> Quantity:
    result_inertia_expr = solve(definition, moment_of_inertia, dict=True)[0][moment_of_inertia]
    result_expr = result_inertia_expr.subs({particle_mass: mass_, spinning_radius: radius_})
    return expr_to_quantity(result_expr, 'moment_of_inertia')