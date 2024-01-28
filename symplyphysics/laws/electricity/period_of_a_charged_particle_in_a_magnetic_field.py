from sympy import (Eq, solve, pi)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Let an arbitrary particle move in a magnetic field around a circle. Then period of its motion depends on
## magnetic induction and on charge and mass of the particle.

## Law is: T = 2 * pi * m / (q * B), where
## T - period,
## m - mass,
## q - charge,
## B - induction.

period = Symbol("period", units.time)

mass = Symbol("mass", units.mass)
charge = Symbol("charge", units.charge)
induction = Symbol("induction", units.magnetic_density)

law = Eq(period, 2 * pi * mass / (charge * induction))


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=mass, charge_=charge, induction_=induction)
@validate_output(period)
def calculate_period(mass_: Quantity, charge_: Quantity, induction_: Quantity) -> Quantity:
    result_period_expr = solve(law, period, dict=True)[0][period]
    result_expr = result_period_expr.subs({mass: mass_, charge: charge_, induction: induction_})
    return Quantity(result_expr)
