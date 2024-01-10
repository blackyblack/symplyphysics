from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The product of the area and the fluid speed, which is called volume flux, is constant 
## at all points along the tube of flow of an incompressible liquid.

# Law: A*v = constant
## A - area of the tube of flow
## v - fluid speed

# It is also called the equation of continuity

# Condition
## The fluid should be ideal, i.e. nonviscous, in laminar flow, incompressible and irrotational.

area_start = Symbol("A_start", units.area)
v_start = Symbol("v_start", units.velocity)
area_end = Symbol("A_end", units.area)
v_end = Symbol("v_end", units.velocity)

law = Eq(area_start * v_start, area_end * v_end)


@validate_input(area_start_=area_start, v_start_=v_start, area_end_=area_end)
@validate_output(v_end)
def calculate_flow_speed(area_start_: Quantity, v_start_: Quantity, area_end_: Quantity) -> Quantity:
    solved = solve(law, v_end, dict=True)[0][v_end]
    value = solved.subs({
        area_start: area_start_,
        v_start: v_start_,
        area_end: area_end_,
    })
    return Quantity(value)
