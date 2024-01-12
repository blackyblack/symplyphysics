from sympy import Eq, dsolve, Derivative
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
## Bernoulli's equation applied to an ideal liquid specifies that the inner 
## pressure of the fluid is constant at all points along a streamline.

# Law: d(P_inner)/dt = 0
## P_inner - inner pressure of fluid at chosen point
## t - time

# Condition
## The fluid must be ideal, i.e.
## 1) nonviscous
## 2) in steady (laminar) flow
## 3) incompressible
## 4) irrotational

time = Symbol("time", units.time)
inner_pressure = Function("inner_pressure", units.pressure)

law = Eq(Derivative(inner_pressure(time), time), 0)


def print_law() -> str:
    return print_expression(law)


@validate_input(inner_pressure_before_=inner_pressure)
@validate_output(inner_pressure)
def calculate_inner_pressure(inner_pressure_before_: Quantity) -> Quantity:
    dsolved = dsolve(law, inner_pressure(time))
    result_expr = dsolved.subs("C1", inner_pressure_before_).rhs
    return Quantity(result_expr)
