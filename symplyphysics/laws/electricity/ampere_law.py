from sympy import (Eq, solve, cos)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, angle_type)

# Description
## Ampere's law is the law that determines the force with which a magnetic field acts on
## a small segment of a conductor with a current. The force turns out to be linearly dependent
## on both current and magnetic induction.

## Law is: F = I * L * B * cos(a), where
## F - force,
## I - current,
## L - length of conductor,
## B - magnetic induction,
## a - angle between magnetic induction and current direction.

# Conditions:
## - Сonductor must be thin;
## - Сonductor must be of finite length.

force = Symbol("force", units.force)

current = Symbol("current", units.current)
length = Symbol("length", units.length)
induction = Symbol("induction", units.magnetic_density)
angle = Symbol("angle", angle_type)

law = Eq(force, current * length * induction * cos(angle))


def print_law() -> str:
    return print_expression(law)


@validate_input(current_=current, length_=length, angle_=angle, induction_=induction)
@validate_output(force)
def calculate_force(current_: Quantity, length_: Quantity, angle_: float | Quantity,
    induction_: Quantity) -> Quantity:
    result_expr = solve(law, force, dict=True)[0][force]
    result_expr = result_expr.subs({
        current: current_,
        length: length_,
        angle: angle_,
        induction: induction_,
    })
    return Quantity(result_expr)
