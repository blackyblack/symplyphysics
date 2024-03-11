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
## When an object undergoes hydraulic compression due to a stress exerted by a surrounding liquid,
## the pressure (hydraulic stress) on the object due to the fluid is proportional to the fractional
## change in the object's volume due to that pressure and the bulk modulus of the object. Thus, bulk
## modulus of a substance is a measure of its resistance to bulk compression.

# Law: p = B * (dV / V)
## p - hydraulic stress (pressure)
## B - bulk modulus of object
## dV - absolute value of change in object's volume
## V - object's initial volume

hydraulic_stress = Symbol("hydraulic_stress", units.pressure)
bulk_modulus = Symbol("bulk_modulus", units.pressure)
change_in_volume = Symbol("change_in_volume", units.volume)
initial_volume = Symbol("initial_volume", units.volume)

law = Eq(hydraulic_stress, bulk_modulus * change_in_volume / initial_volume)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    bulk_modulus_=bulk_modulus,
    change_in_volume_=change_in_volume,
    initial_volume_=initial_volume,
)
@validate_output(hydraulic_stress)
def calculate_hydraulic_stress(
    bulk_modulus_: Quantity,
    change_in_volume_: Quantity,
    initial_volume_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        bulk_modulus: bulk_modulus_,
        change_in_volume: change_in_volume_,
        initial_volume: initial_volume_,
    })
    return Quantity(result)
