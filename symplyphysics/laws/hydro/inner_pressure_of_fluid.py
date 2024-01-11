from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.laws.hydro.dynamic_pressure_from_velocity import dynamic_pressure
from symplyphysics.laws.hydro.hydrostatic_pressure_from_density_and_depth import hydrostatic_pressure

# Description
## Inner pressure of an ideal fluid is the sum of static, dynamic, and hydrostatic pressure at chosen point.

# Law: P_inner = P_static + P_dynamic + P_hydrostatic
## P_inner -- inner pressure
## P_static -- static pressure
## P_dynamic -- dynamic pressure
## P_hydrostatic -- hydrostatic pressure

# Condition: this definition applies to an ideal liquid, namely one that is
## 1) nonviscous
## 2) in steady (laminar) flow
## 3) incompressible
## 4) irrotational

inner_pressure = Symbol("inner_pressure", units.pressure)
static_pressure = Symbol("static_pressure", units.pressure)

law = Eq(inner_pressure, static_pressure + dynamic_pressure + hydrostatic_pressure)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    static_pressure_=static_pressure,
    dynamic_pressure_=dynamic_pressure,
    hydrostatic_pressure_=hydrostatic_pressure,
)
@validate_output(inner_pressure)
def calculate_inner_pressure(
    static_pressure_: Quantity,
    dynamic_pressure_: Quantity,
    hydrostatic_pressure_: Quantity,
) -> Quantity:
    result_expr = solve(law, inner_pressure)[0]
    result_inner_pressure = result_expr.subs({
        static_pressure: static_pressure_,
        dynamic_pressure: dynamic_pressure_,
        hydrostatic_pressure: hydrostatic_pressure_,
    })
    return Quantity(result_inner_pressure)
