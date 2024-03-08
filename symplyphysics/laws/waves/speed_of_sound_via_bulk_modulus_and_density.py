from sympy import Eq, sqrt
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Sound waves are longitudinal mechanical waves that can travel through any type of material,
## be it solid, liquid, or gas. The phase velocity of the a sound wave depends on the bulk modulus 
## of the medium it is traveling in and its density.

# Law: v = sqrt(B / rho)
## v - phase velocity of sound wave
## B - bulk modulus of medium
## rho - density of medium

phase_velocity = Symbol("phase_velocity", units.velocity, positive=True)
bulk_modulus = Symbol("bulk_modulus", units.pressure, real=True)
density = Symbol("density", units.mass / units.volume, positive=True)

law = Eq(phase_velocity, sqrt(bulk_modulus / density))


def print_law() -> str:
    return print_expression(law)


@validate_input(bulk_modulus_=bulk_modulus, density_=density,)
@validate_output(phase_velocity)
def calculate_speed_of_sound(bulk_modulus_: Quantity, density_: Quantity) -> Quantity:
    result = law.rhs.subs({
        bulk_modulus: bulk_modulus_,
        density: density_,
    })
    return Quantity(result)
