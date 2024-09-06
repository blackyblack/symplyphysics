from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The magnetic moment is the main physical quantity characterizing the magnetic properties of a substance,
## that is, the ability to create and perceive a magnetic field.

## Law is: m = I * S, where
## m - magnetic moment,
## I - current,
## S - area of closed contour, which is formed by a closed thin conductor through which current flows.

# Conditions:
## - Ideally, the conductor itself should be infinitely thin. But a conductor can have
##   any shape and size if its dimensions can be neglected relative to the size of the
##   closed circuit that it forms.
## - The plane must be flat for this formula to be applicable. If there are several turns
##   of the conductor, then the formula will be similar, but still different.

moment = Symbol("moment", units.current * units.area)

current = Symbol("current", units.current)
area = Symbol("area", units.area)

law = Eq(moment, current * area)


@validate_input(current_=current, area_=area)
@validate_output(moment)
def calculate_moment(current_: Quantity, area_: Quantity) -> Quantity:
    result_expr = solve(law, moment, dict=True)[0][moment]
    result_expr = result_expr.subs({current: current_, area: area_})
    return Quantity(result_expr)
