from sympy import Eq, Integral
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

# Law: W = Integral(F(x), (x, x_start, x_end))
## W - work done by force F
## F = F(x) - force in question
## x - position of the object on the axis
## x_start, x_end - starting and ending positions of the object during its movement, respectively

work = Symbol("work", units.energy)
force = Function("force", units.force)
position = Symbol("position", units.length)
position_start = Symbol("position_start", units.length)
position_end = Symbol("position_end", units.length)

law = Eq(work, Integral(force(position), (position, position_start, position_end)))


def print_law() -> str:
    return print_expression(law)


# Assuming the force changes linearly with respect to position
@validate_input(
    force_start_=force,
    force_end_=force,
    position_start_=position,
    position_end_=position,
)
@validate_output(work)
def calculate_work(
    force_start_: Quantity,
    force_end_: Quantity,
    position_start_: Quantity,
    position_end_: Quantity,
) -> Quantity:
    # Using the two-point line equation: (y - y1)/(y2 - y1) = (x - x1)/(x2 - x1)
    force_function = (
        (force_end_ - force_start_)
        * (position - position_start_)
        / (position_end - position_start_)
        + force_start_
    )
    result = law.rhs.subs({
        force(position): force_function,
        position_start: position_start_,
        position_end: position_end_,
    })
    result_work = result.doit()
    return Quantity(result_work)
