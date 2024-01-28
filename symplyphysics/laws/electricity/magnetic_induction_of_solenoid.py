from sympy import (Eq, solve)
from sympy.physics.units import magnetic_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless
)

# Description
## A solenoid is a cylindrical coil consisting of a large number of turns of wire forming a helical line.
## The magnetic induction of solenoid depends on current, number of turns, length of solenoid and material.

## Law is: B = mu * mu0 * I * N / l, where
## B - induction,
## mu - relative permeability of the core inside of a solenoid,
## mu0 - magnetic constant,
## I - current,
## N - number of turns,
## l - length of solenoid.

induction = Symbol("induction", units.magnetic_density)

current = Symbol("current", units.current)
length = Symbol("length", units.length)
number_turns = Symbol("number_turns", dimensionless)
relative_permeability = Symbol("relative_permeability", dimensionless)

law = Eq(induction, relative_permeability * magnetic_constant * current * number_turns / length)


def print_law() -> str:
    return print_expression(law)


@validate_input(current_=current,
    length_=length,
    relative_permeability_=relative_permeability,
    number_turns_=number_turns)
@validate_output(induction)
def calculate_induction(current_: Quantity, length_: Quantity, relative_permeability_: float, number_turns_: float) -> Quantity:
    #TODO: throw on negative number of turns
    result_expr = solve(law, induction, dict=True)[0][induction]
    result_expr = result_expr.subs({
        current: current_,
        length: length_,
        relative_permeability: relative_permeability_,
        number_turns: number_turns_,
    })
    return Quantity(result_expr)
