from sympy import Eq, solve, dsolve, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Bernoulli's equation applied to an ideal liquid specifies that the sum of static,
## dynamic and hydrostatic pressure is constant at all points along a streamline.

# Law: d[P(t)/rho + v(t)**2/2 + g*h(t)]/dt = 0
## P - pressure at chosen point
## rho - density of the fluid
## v - flow speed at chosen point
## g - gravitational constant
## h - elevation of chosen point above reference plane

# Condition
## The fluid must be ideal, i.e.
## 1) nonviscous
## 2) in steady (laminar) flow
## 3) incompressible
## 4) irrotational

time = Symbol("t", units.time)
pressure = Function("P", units.pressure)
density = Symbol("rho", units.mass / units.volume)
speed = Function("v", units.speed)
elevation = Function("h", units.length)

law = Eq(Derivative(
    pressure(time) / density 
    + speed(time) ** 2 / 2 
    + units.acceleration_due_to_gravity * elevation(time)
, time), 0)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    pressure_before_=pressure, density_=density, speed_before_=speed,
    elevation_before_=elevation, speed_after_=speed, elevation_after_=elevation,
)
@validate_output(pressure)
def calculate_pressure(
    pressure_before_: Quantity, density_: Quantity, speed_before_: Quantity,
    elevation_before_: Quantity, speed_after_: Quantity, elevation_after_: Quantity,
) -> Quantity:
    dsolved = dsolve(law, pressure(time))
    c1_value = solve(dsolved, "C1")[0].subs({
        pressure(time): pressure_before_,
        density: density_,
        speed(time): speed_before_,
        elevation(time): elevation_before_,
    })
    result_expr = dsolved.subs({
        "C1": c1_value,
        density: density_,
        speed(time): speed_after_,
        elevation(time): elevation_after_,
    }).rhs
    return Quantity(result_expr)
