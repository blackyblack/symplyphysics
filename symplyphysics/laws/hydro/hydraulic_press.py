from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output)

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


first_force = Symbol("first_force", units.force)
first_area = Symbol("first_area", units.area)
second_force = Symbol("second_force", units.force)
second_area = Symbol("second_area", units.area)

law = Eq(first_force / first_area, second_force / second_area)


def print_law() -> str:
    return print_expression(law)


@validate_input(first_force_=first_force, first_area_=first_area, second_area_=second_area)
@validate_output(second_force)
def calculate_force(first_force_: Quantity, first_area_, second_area_: Quantity) -> Quantity:
    result_expr = solve(law, second_force, dict=True)[0][second_force]
    result_force = result_expr.subs({
        first_force: first_force_,
        first_area: first_area_,
        second_area: second_area_,
    })

    return Quantity(result_force)