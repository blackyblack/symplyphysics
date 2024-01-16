from sympy import Eq, solve, Abs
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The magnitude of the electric field set up by a point charge is linearly proportional 
## to the absolute value of the charge and the square inverse of the distance to it.

# Law: E = k_e * |q| / r^2
## E - magnitude of electric field
## k_e - Coulomb's constant
## q - charge of point charge
## r - distance to point charge

electric_field = Symbol("electric_field", units.force / units.charge)
point_charge = Symbol("point_charge", units.charge)
distance = Symbol("distance", units.length)

law = Eq(electric_field, units.coulomb_constant * Abs(point_charge) / distance**2)


def print_law() -> str:
    return print_expression(law)


@validate_input(point_charge_=point_charge, distance_=distance)
@validate_output(electric_field)
def calculate_electric_field(point_charge_: Quantity, distance_: Quantity) -> Quantity:
    result = solve(law, electric_field)[0]
    result_field = result.subs({
        point_charge: point_charge_,
        distance: distance_,
    })
    return Quantity(result_field)
