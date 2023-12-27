from sympy import sqrt
from sympy import Eq, solve
from sympy.physics.units import speed_of_light, convert_to
from symplyphysics import units, Quantity, Symbol, print_expression, validate_input, validate_output

# Description
## The relativistic energy of a particle of rest mass m moving in 
# your frame of reference at speed v 

# Law: E = gamma * m * speed_of_light**2, where
## E is energy of the particle
## m is the rest mass of the particle
## speed_of_light is the speed of light.
## v is the speed of the particle
## gamma is 1 / sqrt(1 - (v/speed_of_light) ** 2)

particle_energy = Symbol("particle_energy", units.energy)
rest_mass = Symbol("rest_mass", units.mass)
particle_velocity = Symbol("particle_velocity", units.velocity)
law = Eq(particle_energy, 1 / sqrt(1 - (particle_velocity / speed_of_light) ** 2) * rest_mass * speed_of_light ** 2)


def print_law() -> str:
    return print_expression(law)


@validate_input(rest_mass_=rest_mass, particle_veclocity_=particle_velocity)
@validate_output(particle_energy)
def calculate_particle_energy(rest_mass_: Quantity,
                              particle_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, particle_energy, dict=True)[0][particle_energy]
    particle_energy_res = result_expr.subs({rest_mass: rest_mass_, particle_velocity: particle_velocity_})
    return Quantity(particle_energy_res)
