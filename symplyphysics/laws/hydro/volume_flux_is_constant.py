from sympy import Eq, dsolve, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)

# Description
## The product of the area and the fluid speed, which is called volume flux, is constant 
## at all points along the tube of flow of an incompressible liquid.

# Law: d(A*v)/dt = 0
## A - area of the tube of flow
## v - fluid speed

# It is also called the equation of continuity

# Condition
## The fluid should be ideal, i.e.
## 1) nonviscous,
## 2) in steady (laminar) flow,
## 3) incompressible,
## 4) irrotational.

time = Symbol("time", units.time)
tube_area = Function("tube_area", units.area)
fluid_speed = Function("fluid_speed", units.velocity)

law = Eq(Derivative(tube_area(time) * fluid_speed(time), time), 0)


@validate_input(tube_area_before_=tube_area, fluid_speed_before_=fluid_speed, tube_area_after_=tube_area)
@validate_output(units.velocity)
def calculate_flow_speed(tube_area_before_: Quantity, fluid_speed_before_: Quantity, tube_area_after_: Quantity) -> Quantity:
    solved = dsolve(law, fluid_speed(time))
    result_expr = solved.subs({
        "C1": tube_area_before_ * fluid_speed_before_,
        tube_area(time): tube_area_after_
    }).rhs
    return Quantity(result_expr)
