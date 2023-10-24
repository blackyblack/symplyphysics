from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

# Description
## Power factor is property of any AC consumer. Commonly not all power consumed from source makes useful work.
## The part of consumed power which makes work called active power.
## Part of it may be returned to source and it is also known as reactive power.
## Law: Q = P / S, where
## Q is power factor,
## P is active power consumed,
## S is total power consumed.

full_power = Symbol("full_power", units.power)
active_power = Symbol("circular_frequency", units.power)
power_factor = Symbol("power_factor", dimensionless)

law = Eq(power_factor, active_power / full_power)

def print_law() -> str:
    return print_expression(law)


@validate_input(active_power_=active_power, full_power_=full_power)
@validate_output(power_factor)
def calculate_power_factor(active_power_: Quantity, full_power_: Quantity) -> Quantity:
    result_factor_expr = solve(law, power_factor, dict=True)[0][power_factor]
    result_expr = result_factor_expr.subs({
        active_power: active_power_,
        full_power: full_power_
    })
    return Quantity(result_expr)
