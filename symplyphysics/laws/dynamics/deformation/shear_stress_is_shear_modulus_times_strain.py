from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## When an object is under a shearing stress, the shear stress applied to the object is proportional
## to the shear modulus of the object and the ratio between the displacement of one end of the object
## in the direction of the applied (shear) force and the length of the object.

# Law: tau = G * (dx / L)
## tau - [shear stress](../../hydro/shear_stress_is_shear_force_over_area.py)
## G - shear modulus of the object
## dx - transverse displacement, i.e. displacement of one end of the object in the direction of the force
## L - initial length of the area which the force is applied to

shear_stress = Symbol("shear_stress", units.pressure)
shear_modulus = Symbol("shear_modulus", units.pressure)
transverse_displacement = Symbol("transverse_displacement", units.length)
initial_length = Symbol("initial_length", units.length)

law = Eq(shear_stress, shear_modulus * transverse_displacement / initial_length)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    shear_modulus_=shear_modulus,
    transverse_displacement_=transverse_displacement,
    initial_length_=initial_length,
)
@validate_output(shear_stress)
def calculate_shear_stress(
    shear_modulus_: Quantity,
    transverse_displacement_: Quantity,
    initial_length_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        shear_modulus: shear_modulus_,
        transverse_displacement: transverse_displacement_,
        initial_length: initial_length_,
    })
    return Quantity(result)
