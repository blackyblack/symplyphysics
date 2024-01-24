from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output)

# Description
## Linear charge density is the quantity of charge per unit length, at any point on a line charge distribution.
## Charge density can be either positive or negative, since electric charge can be either positive or negative.


## Law: λ = q / l
## Where:
## λ is linear charge density
## q is charge
## l is the length over which charge is distributed


linear_charge_density = Symbol("linear_charge_density", units.charge / units.length)
charge = Symbol("charge", units.charge)
length = Symbol("length", units.length)

law = Eq(linear_charge_density, charge / length)


def print_law() -> str:
    return print_expression(law)


@validate_input(charge_=charge, length_=length)
@validate_output(linear_charge_density)
def calculate_linear_charge_density(charge_, length_: Quantity) -> Quantity:
    result_expr = solve(law, linear_charge_density, dict=True)[0][linear_charge_density]
    result_linear_charge_density = result_expr.subs({
        charge: charge_,
        length: length_,
    })
    return Quantity(result_linear_charge_density)
