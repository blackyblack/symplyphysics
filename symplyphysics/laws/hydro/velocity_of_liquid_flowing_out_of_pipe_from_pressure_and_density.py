from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.hydro import hydrostatic_pressure_from_density_and_depth as pressure_law
from symplyphysics.laws.hydro import velocity_from_height as velocity_law

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

law = Eq(velocity, sqrt(2 * pressure / density))

# This law might be derived via "hydrostatic_pressure_from_density_and_depth" law
# and "velocity_from_height" law.

pressure_law_applied = pressure_law.law.subs({
    pressure_law.density: density,
    pressure_law.hydrostatic_pressure: pressure
})
height_derived = solve(pressure_law_applied, pressure_law.depth, dict=True)[0][pressure_law.depth]

velocity_law_applied = velocity_law.law.subs({
    velocity_law.height_above_hole: height_derived
})
velocity_derived = solve(velocity_law_applied, velocity_law.liquid_velocity, dict=True)[0][velocity_law.liquid_velocity]

# Check if derived velocity is same as declared.
assert expr_equals(velocity_derived, law.rhs)


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
