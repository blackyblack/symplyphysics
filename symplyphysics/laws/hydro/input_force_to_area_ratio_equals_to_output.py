from sympy import (Eq, solve, dsolve)
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol, print_expression,
    validate_input, validate_output)
from symplyphysics.laws.dynamics import pressure_from_force_and_area as pressure_law
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.hydro import inner_pressure_is_constant as constant_pressure_law

# Description
## If both vertically positioned cylinders of communicating vessels are closed with pistons, then with the help of external forces applied to the pistons,
## a large pressure can be created in the liquid, many times exceeding the hydrostatic pressure at any point in the system.
## If the pistons have different areas, then different forces act on them from the liquid side. The same modulus, but oppositely directed external forces must be applied to the pistons to keep the system in balance.

## Law: F1 / S1 = F2 / S2
## Where:
## F1 is the force acting on the first piston
## F2 is the force acting on the second piston
## S1 is the area of the first piston
## S2 is the area of the second piston

## Conditions
## This ratio is performed only in an ideal hydraulic press, i.e. one in which there is no friction.

input_force = clone_symbol(symbols.dynamics.force)
input_area = Symbol("input_area", units.area)
output_force = clone_symbol(symbols.dynamics.force)
output_forces_area = Symbol("output_forces_area", units.area)

law = Eq(input_force / input_area, output_force / output_forces_area)

pressure_input = pressure_law.law.rhs.subs({
    pressure_law.force: input_force,
    pressure_law.area: input_area
})

pressure_output = pressure_law.law.rhs.subs({
    pressure_law.force: output_force,
    pressure_law.area: output_forces_area,
})

## If the pistons are in equilibrium, then the pressures pressure_input and pressure_output are equal
dsolved = dsolve(constant_pressure_law.law,
    constant_pressure_law.inner_pressure(constant_pressure_law.time))
dsolved_input = dsolved.subs(constant_pressure_law.inner_pressure(constant_pressure_law.time),
    pressure_input)
dsolved_output = dsolved.subs(constant_pressure_law.inner_pressure(constant_pressure_law.time),
    pressure_output)
solved_input = solve([dsolved_input, dsolved_output], (pressure_input, "C1"),
    dict=True)[0][pressure_input]
pressure_equation = Eq(pressure_input, solved_input)

assert expr_equals(law.rhs, pressure_equation.rhs)
assert expr_equals(law.lhs, pressure_equation.lhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(input_force_=input_force,
    input_area_=input_area,
    output_forces_area_=output_forces_area)
@validate_output(output_force)
def calculate_output_force(input_force_: Quantity, input_area_: Quantity,
    output_forces_area_: Quantity) -> Quantity:
    result_expr = solve(law, output_force, dict=True)[0][output_force]
    result_force = result_expr.subs({
        input_force: input_force_,
        input_area: input_area_,
        output_forces_area: output_forces_area_,
    })
    return Quantity(result_force)
