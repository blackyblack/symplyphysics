from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Consider a tube through which a liquid flows under arbitrary pressure. Then the velocity of the
## liquid flowing out of the pipe will depend only on the pressure and density of the liquid.

## Law is: v = sqrt(2 * p / p0), where
## v - velocity of the liquid flowing out of the pipe,
## p - pressure,
## p0 - density of the liquid.

velocity = Symbol("velocity", units.velocity)

pressure = Symbol("pressure", units.pressure)
density = Symbol("density", units.mass / units.volume)

law = Eq(velocity, sqrt( 2 * pressure / density))


def print_law() -> str:
    return print_expression(law)


@validate_input(pressure_=pressure, density_=density)
@validate_output(velocity)
def calculate_velocity(pressure_: Quantity, density_: Quantity) -> Quantity:
    result_expr = solve(law, velocity, dict=True)[0][velocity]
    result_expr = result_expr.subs({
        pressure: pressure_,
        density: density_,
    })
    return Quantity(result_expr)
