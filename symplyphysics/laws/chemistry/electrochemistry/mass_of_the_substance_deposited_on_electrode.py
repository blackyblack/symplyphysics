from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

# Description
## Faraday's first law of electrolysis: the mass of a substance deposited on an electrode during electrolysis
## is directly proportional to the amount of electricity transferred to this electrode. By the amount of
## electricity, we mean the total electric charge that has passed through the surface of the electrode.

## Law is: m = k * I * t, where
## m - mass,
## k - electrochemical equivalent,
## I - current,
## t - time.

equivalent = Symbol("equivalent", units.mass / units.charge)
current = Symbol("current", units.current)
time = Symbol("time", units.time)
mass_deposited = symbols.mass

law = Eq(mass_deposited, equivalent * current * time)


def print_law() -> str:
    return print_expression(law)


@validate_input(equivalent_=equivalent, current_=current, time_=time)
@validate_output(mass_deposited)
def calculate_mass(equivalent_: Quantity, current_: Quantity, time_: Quantity) -> Quantity:
    result_expr = solve(law, mass_deposited, dict=True)[0][mass_deposited]
    result_expr = result_expr.subs({equivalent: equivalent_, current: current_, time: time_})
    return Quantity(result_expr)
