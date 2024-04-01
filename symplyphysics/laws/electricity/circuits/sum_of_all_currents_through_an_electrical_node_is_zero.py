from typing import Sequence
from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, print_expression, validate_input, validate_output,
    SymbolIndexed, SumIndexed, global_index)

# Description
## sum(I) = 0
# Where:
## I is a current flowing through electrical node.

## In other words, as electrical charge is being neither created nor accumulated in electrical node,
## sum of all input currents are equal to sum of output current. If we assert input current is positive and output is negative, we gain summary current as 0.
## This property of electrical node is called Kirchhoff law #1.
## Assert there are minimum 2 currents flowing through any node

current = SymbolIndexed("current", units.current)
law = Eq(SumIndexed(current[global_index], global_index), 0)


def print_law() -> str:
    return print_expression(law)


@validate_input(currents_=current)
@validate_output(units.current)
def calculate_current_from_array(currents_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(currents_) + 1))
    currents_law = law.subs(global_index, local_index)
    currents_law = currents_law.doit()
    unknown_current = current[len(currents_) + 1]
    solved = solve(currents_law, unknown_current, dict=True)[0][unknown_current]
    for i, v in enumerate(currents_):
        solved = solved.subs(current[i + 1], v)
    return Quantity(solved)
