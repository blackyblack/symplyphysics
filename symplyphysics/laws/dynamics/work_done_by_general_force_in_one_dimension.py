from sympy import Eq, Integral, symbols, solve
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
## Assuming one-dimensional environment, when the force F on a particle-like object depends
## on the position of the object, the work done by F on the object while the object moves
## from one position to another is to be found by integrating the force along the path of the
## object.

# Law: W = Integral(F(x), x)
## W - work done by force F
## F = F(x) - force in question
## x - position of the object on the axis

work = Symbol("work", units.energy)
force = Function("force", units.force)
position = Symbol("position", units.length)

law = Eq(work, Integral(force(position), position))


def print_law() -> str:
    return print_expression(law)


# Assuming the force changes linearly with respect to position
@validate_input(
    force_before_=force,
    force_after_=force,
    position_before_=position,
    position_after_=position,
)
@validate_output(work)
def calculate_work(
    force_before_: Quantity,
    force_after_: Quantity,
    position_before_: Quantity,
    position_after_: Quantity,
) -> Quantity:
    a, b = symbols("a b")
    solved = solve(
        [
            Eq(force_before_, a * position_before_ + b),
            Eq(force_after_, a * position_after_ + b),
        ],
        (a, b),
        dict=True,
    )[0]
    force_function = solved[a] * position + solved[b]
    expr = law.rhs.subs(force(position), force_function).doit()
    value = expr.subs(position, position_after_) - expr.subs(position, position_before_)
    return Quantity(value)
