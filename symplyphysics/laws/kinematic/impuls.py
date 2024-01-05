from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression,
                           validate_input, validate_output)


# Descroption
## The momentum of a body (the amount of motion) is a vector physical quantity, the modulus of which is equal to the product of the mass of the body by its speed.
## p = m * v
## Where
## m - mass of the body
## v - velocity

velocity = Symbol("velocity", units.velocity)
mass_of_the_body = Symbol("mass_of_the_body", units.mass)

impuls = Symbol("impuls", units.velocity * units.mass)

law = Eq(impuls, mass_of_the_body * velocity)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_of_the_body_=mass_of_the_body, velocity_=velocity)
@validate_output(velocity)
def calculate_velocity(mass_of_the_body_: Quantity, velocity_: Quantity) -> Quantity:
    result_impuls = solve(law, velocity, dict=True)[0][velocity]
    result_expr = result_impuls.subs({mass_of_the_body: mass_of_the_body_, velocity: velocity_})
    return Quantity(result_expr)
