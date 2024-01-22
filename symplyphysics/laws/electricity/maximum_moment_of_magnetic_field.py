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
## The magnetic moment is the main physical quantity characterizing the magnetic properties of a substance,
## that is, the ability to create and perceive a magnetic field.

## Law is: m = I * S, where
## m - magnetic moment,
## I - current,
## S - area of the circuit through which the current flows.

moment = Symbol("moment", units.current * units.area)

current = Symbol("current", units.current)
area = Symbol("area", units.area)

law = Eq(moment, current * area)


def print_law() -> str:
    return print_expression(law)


@validate_input(current_=current, area_=area)
@validate_output(moment)
def calculate_moment(current_: Quantity, area_: Quantity) -> Quantity:
    result_expr = solve(law, moment, dict=True)[0][moment]
    result_expr = result_expr.subs({
        current: current_,
        area: area_
    })
    return Quantity(result_expr)
