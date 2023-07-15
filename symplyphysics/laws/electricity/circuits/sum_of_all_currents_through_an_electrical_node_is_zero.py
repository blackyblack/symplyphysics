from sympy import (Eq, solve)
from symplyphysics import (Symbol, units, Quantity, print_expression, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

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


def print_law() -> str:
    return print_expression(law)


@validate_input(currents_=currents)
@validate_output(units.current)
def calculate_current_from_array(currents_: list[Quantity]) -> Quantity:
    current_symbols = tuple_of_symbols("current", units.current, len(currents_) + 1)
    unknown_current = current_symbols[len(currents_)]
    currents_law = law.subs(currents, current_symbols).doit()
    solved = solve(currents_law, unknown_current, dict=True)[0][unknown_current]
    for (from_, to_) in zip(current_symbols, currents_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
