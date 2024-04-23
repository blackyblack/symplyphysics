from sympy import (Eq, solve)
from symplyphysics import (symbols, units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Pressure has a direct relationship with force. Assuming that the area is constant, pressure increases as the force applied also increases.

## Law: P = F / A
## Where:
## P is pressure
## F is the force
## A is the area

pressure = Symbol("pressure", units.pressure)
force = symbols.dynamics.force
area = Symbol("area", units.area)

law = Eq(pressure, force / area)


def print_law() -> str:
    return print_expression(law)


@validate_input(force_=force, area_=area)
@validate_output(pressure)
def calculate_pressure(force_: Quantity, area_: Quantity) -> Quantity:
    result_expr = solve(law, pressure, dict=True)[0][pressure]
    result_pressure = result_expr.subs({
        force: force_,
        area: area_,
    })

    return Quantity(result_pressure)
