from sympy import (Eq, solve)
from symplyphysics import (
    units, Quantity, Symbol, print_expression, validate_input, validate_output, dimensionless)


diameter = Symbol("d", units.length)  # hydraulic diameter of the pipe
# density of the fluid kg/m^3
density = Symbol("rho", units.mass / units.volume)
velocity = Symbol("v", units.length / units.time)  # velocity of the fluid m/s
# dynamic viscosity of the fluid Pa*s
dynamic_viscosity = Symbol("mu", units.pressure * units.time)
reynolds_number = Symbol("Re", dimensionless)

law = Eq(reynolds_number, (density * velocity * diameter) / dynamic_viscosity)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    diameter_=diameter,
    density_=density,
    velocity_=velocity,
    dynamic_viscosity_=dynamic_viscosity
)
@validate_output(reynolds_number)
def calculate_reynolds_number(diameter_, density_, velocity_, dynamic_viscosity_):
    result_reynolds_number_expr = solve(law, reynolds_number, dict=True)[
        0][reynolds_number]
    result_expr = result_reynolds_number_expr.subs({
        diameter: diameter_,
        density: density_,
        velocity: velocity_,
        dynamic_viscosity: dynamic_viscosity_
    })
    return Quantity(result_expr, dimension=dimensionless)
