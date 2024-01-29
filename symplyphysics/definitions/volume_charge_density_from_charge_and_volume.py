from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output)
from symplyphysics.definitions import volume_number_density
from symplyphysics.laws.electricity import charge_is_quantized

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

definition = Eq(volume_charge_density, charge / volume)

definition_units_SI = units.coulomb / units.meter**3

## Proof

## Let's express the number of objects from the volume number density definition
## In this case, the objects are charged particles.
expr = volume_number_density.definition.subs({
        volume_number_density.volume: volume,
        volume_number_density.objects: charge
    })

derived_law = solve(expr,  volume_number_density.number_density, dict=True)[0][volume_number_density.number_density]
assert definition.rhs == derived_law

def print_definition() -> str:
    return print_expression(definition)


@validate_input(charge_=charge, volume_=volume)
@validate_output(volume_charge_density)
def calculate_volume_charge_density(charge_: Quantity, volume_: Quantity) -> Quantity:
    result_expr = solve(definition, volume_charge_density, dict=True)[0][volume_charge_density]
    result_volume_charge_density = result_expr.subs({
        charge: charge_,
        volume: volume_,
    })
    return Quantity(result_volume_charge_density)
