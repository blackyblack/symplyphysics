from sympy import (Eq, solve, pi)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Brownian motion is the random motion of microscopic visible suspended particles of a solid substance in a liquid or gas caused by the thermal motion of particles of a liquid or gas.
## The diffusion coefficient is a quantitative characteristic of the diffusion rate, equal to the amount of matter passing per unit time through a section of a unit area as a result of the thermal motion of molecules with a concentration gradient equal to one(corresponding to a change of 1 mol/l → 0 mol/l per unit length).
## The diffusion coefficient is determined by the properties of the medium and the type of diffusing particles.
## This law is also known as Stokes–Einstein–Sutherland relation.

## Law: D = R * T / (6 * Na * pi * r * eta)
## Where:
## D is diffusion coefficient
## R is universal gas constant
## T is temperature
## Na is avogadro number
## r is particle radius
## eta is dynamic viscosity

# Conditions:
## - Particle displacements in any direction are equally likely
## - The inertia of a Brownian particle can be neglected compared to the influence of friction forces
## - Particles are spherical
## - Low Reynolds number, ie non turbulent flow

diffusion_coefficient = Symbol("diffusion_coefficient", units.area / units.time)
temperature = Symbol("temperature", units.temperature)
particle_radius = Symbol("particle_radius", units.length)
dynamic_viscosity = Symbol("dynamic_viscosity", units.pressure * units.time)

law = Eq(
    diffusion_coefficient, units.molar_gas_constant * temperature /
    (6 * units.avogadro * pi * particle_radius * dynamic_viscosity))


def print_law() -> str:
    return print_expression(law)


@validate_input(temperature_=temperature,
    particle_radius_=particle_radius,
    dynamic_viscosity_=dynamic_viscosity)
@validate_output(diffusion_coefficient)
def calculate_diffusion_coefficient(temperature_: Quantity, particle_radius_: Quantity,
    dynamic_viscosity_: Quantity) -> Quantity:
    result_expr = solve(law, diffusion_coefficient, dict=True)[0][diffusion_coefficient]
    result_diffusion_coefficient = result_expr.subs({
        temperature: temperature_,
        particle_radius: particle_radius_,
        dynamic_viscosity: dynamic_viscosity_
    })
    return Quantity(result_diffusion_coefficient)
