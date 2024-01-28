from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Surface charge density is the amount of charge per unit of a two-dimensional surface area. It is a measure of how much quantity of electric charge is accumulated over a surface.
## Charge density can be either positive or negative, since electric charge can be either positive or negative.

## Definition: σ = q / s
## Where:
## σ is surface charge density
## q is charge
## s is the area over which charge is distributed

surface_charge_density = Symbol("surface_charge_density", units.charge / units.area)
charge = Symbol("charge", units.charge)
area = Symbol("area", units.area)

definition = Eq(surface_charge_density, charge / area)

definition_units_SI = units.coulomb / units.meter**2


def print_definition() -> str:
    return print_expression(definition)


@validate_input(charge_=charge, area_=area)
@validate_output(surface_charge_density)
def calculate_surface_charge_density(charge_: Quantity, area_: Quantity) -> Quantity:
    result_expr = solve(definition, surface_charge_density, dict=True)[0][surface_charge_density]
    result_surface_charge_density = result_expr.subs({
        charge: charge_,
        area: area_,
    })
    return Quantity(result_surface_charge_density)
