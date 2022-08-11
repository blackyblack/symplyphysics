from sympy import Idx, IndexedBase, Sum
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## sum(I) = 0
## Where I is a current flowing through electrical node.
## In other words, as electrical charge is being neither created nor accumulated in electrical node,
## sum of all input currents are equal to sum of output current. If we asset input current is positive and output is negative, we gain summary current as 0.
## This property of electrical node is called Kirchhoff law #1.
## Assert there are minimum 2 currents flowing through any node

current = IndexedBase('current')
currents_total = symbols('currents_total')
i = symbols('i', cls=Idx)

law = Eq(Sum(current[i], (i, 1, currents_total)), 0)

def print():
    return pretty(law, use_unicode=False)

@validate_input(current_in = units.current)
@validate_output(units.current)
def calculate_current_1(current_in: Quantity) -> Quantity:
    two_currents_law = law.subs(currents_total, 2).doit()
    solved = solve(two_currents_law, current[1], dict=True)[0][current[2]]
    result_expr = solved.subs(current[2], current_in)
    return expr_to_quantity(result_expr, 'current_out')

