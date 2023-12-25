from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input, validate_output)

# Description
## Every liquid has pressure due to its own weight, also known as static pressure. Static and dynamic pressure added together form total pressure of the liquid.
## Law: p = ro * g * h, where
## ro is density of liquid,
## g is gravity at the surface of liquid
## h is height of liquid column or depth within a substance

# Conditions
## Constant density throughout the liquid
## No variation of g, because height of liquid column is often reasonably small compared to the radius of the Earth

liquid_density = Symbol("liquid_density", units.mass / units.volume)
gravitational_acceleration = Symbol("gravity", units.meter / units.second**2) #normal gravity at the equator approx. = 9.7803267715 m/s^2 
height_of_column = Symbol("height", units.meter)
static_pressure = Symbol("static_pressure", units.pressure)

law = Eq(static_pressure, liquid_density * gravity * height)


def print_law() -> str:
    return print_expression(law)


@validate_input(density_=liquid_density, acceleration_=gravitational_acceleration, height_=height_of_column)
@validate_output(static_pressure)
def calculate_pressure(density_: Quantity, height_: Quantity, acceleration_: Quantity) -> Quantity:
    result_pressure_expr = solve(law, static_pressure, dict=True)[0][static_pressure]
    result_expr = result_pressure_expr.subs({liquid_density: density_, gravitational_acceleration: acceleration_, height_of_column: height_})
    return Quantity(result_expr)
