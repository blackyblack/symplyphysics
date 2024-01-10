from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Current density is the amount of charge per unit time that flows through a unit area of a chosen
## cross section. The current density vector is defined as a vector whose magnitude is the electric
## current per cross-sectional area at a given point in space, its direction being that of the motion
## of the positive charges at this point. In SI base units, the electric current density is measured
## in amperes per square metre.

## Law is: j = q * n * v, where
## j - current density of charge carriers,
## q - charge of the charge carrier (electrons or holes),
## n - concentration of charge carriers,
## v - drift velocity of charge carriers (electrons or holes).

density_current = Symbol("density_current", units.current / units.area)

charge_carriers_concentration = Symbol("charge_carriers_concentration", 1 / units.volume)
drift_velocity = Symbol("drift_velocity", units.velocity)
charge = Symbol("charge", units.charge)

law = Eq(density_current, charge * charge_carriers_concentration * drift_velocity)


def print_law() -> str:
    return print_expression(law)


@validate_input(charge_carriers_concentration_=charge_carriers_concentration, drift_velocity_=drift_velocity, charge_=charge)
@validate_output(density_current)
def calculate_current(charge_carriers_concentration_: Quantity, drift_velocity_: Quantity, charge_: Quantity) -> Quantity:
    result_expr = solve(law, density_current, dict=True)[0][density_current]
    result_expr = result_expr.subs({
        charge_carriers_concentration: charge_carriers_concentration_,
        drift_velocity: drift_velocity_,
        charge: charge_
    })
    return Quantity(result_expr)
