from typing import List
from sympy import (Eq, solve)
from symplyphysics import (Symbol, units, expr_to_quantity, Quantity, print_expression, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.quantities import Quantity

# Description
## sum(I) = 0
# Where:
## I is a current flowing through electrical node.

## In other words, as electrical charge is being neither created nor accumulated in electrical node,
## sum of all input currents are equal to sum of output current. If we assert input current is positive and output is negative, we gain summary current as 0.
## This property of electrical node is called Kirchhoff law #1.
## Assert there are minimum 2 currents flowing through any node

currents = Symbol("currents", units.current)
law = Eq(SumArray(currents), 0, evaluate=False)


def print() -> str:
    return print_expression(law)


@validate_input(currents_=currents)
@validate_output(units.current)
def calculate_current_from_array(currents_: List[Quantity]) -> Quantity:
    current_symbols = tuple(Symbol("current" + str(i), units.current) for i in range(len(currents_) + 1))
    unknown_current = current_symbols[len(currents_)]
    currents_law = law.subs(currents, tuple(current_symbols)).doit()
    solved = solve(currents_law, unknown_current, dict=True)[0][unknown_current]
    for idx in range(len(current_symbols) - 1):
        solved = solved.subs(current_symbols[idx], currents_[idx])
    return expr_to_quantity(solved)
