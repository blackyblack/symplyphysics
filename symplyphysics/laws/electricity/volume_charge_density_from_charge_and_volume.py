from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.definitions import number_density_is_number_of_objects_per_unit_volume

# Description
## The volume charge density is defined as the amount of charge present over a unit volume of the conductor.
## Charge density can be either positive or negative, since electric charge can be either positive or negative.

## Definition: ρ = q / V
## Where:
## ρ is volume charge density
## q is charge
## V is the volume over which charge is distributed

volume_charge_density = Symbol("volume_charge_density", units.charge / units.volume)
charge = Symbol("charge", units.charge)
volume = Symbol("volume", units.volume)

law = Eq(volume_charge_density, charge / volume)

## Proof

## Let's express the number of objects from the volume number density definition
## In this case, the objects are charged particles.
expr = number_density_is_number_of_objects_per_unit_volume.definition.subs({
    number_density_is_number_of_objects_per_unit_volume.volume: volume,
    number_density_is_number_of_objects_per_unit_volume.number_of_objects: charge
})

derived_law = solve(expr, number_density_is_number_of_objects_per_unit_volume.number_density,
    dict=True)[0][number_density_is_number_of_objects_per_unit_volume.number_density]
assert law.rhs == derived_law


def print_definition() -> str:
    return print_expression(law)


@validate_input(charge_=charge, volume_=volume)
@validate_output(volume_charge_density)
def calculate_volume_charge_density(charge_: Quantity, volume_: Quantity) -> Quantity:
    result_expr = solve(law, volume_charge_density, dict=True)[0][volume_charge_density]
    result_volume_charge_density = result_expr.subs({
        charge: charge_,
        volume: volume_,
    })
    return Quantity(result_volume_charge_density)
