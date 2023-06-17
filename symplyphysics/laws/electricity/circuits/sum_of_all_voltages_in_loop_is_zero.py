from typing import List
from sympy import (Eq, solve, symbols, Idx, IndexedBase, Sum)
from symplyphysics import (units, expr_to_quantity, Quantity, print_expression,
    validate_input, validate_output)

# Description
## sum(U) = 0
## Where U is a voltage on the element of an electrical loop.
## Loop is a closed sub-circuit in an electrical circuit.
## In other words, sum of all voltage sources in the loop equals to sum of all voltage consumers in this loop.
## This property of electrical loop is called Kirchhoff law #2.
## This law also demonstrates that work to move the unit of charge in electrical field along the closed path is zero.

voltage = IndexedBase("voltage")
voltages_total = symbols("voltages_total")
i = symbols("i", cls=Idx)

law = Eq(Sum(voltage[i], (i, 1, voltages_total)), 0)


def print() -> str:
    return print_expression(law)


@validate_input(voltages_=units.voltage)
@validate_output(units.voltage)
def calculate_voltage(voltages_: List[Quantity]) -> Quantity:
    voltages_law = law.subs(voltages_total, len(voltages_) + 1).doit()
    unknown_voltage = voltage[len(voltages_) + 1]
    solved = solve(voltages_law, unknown_voltage, dict=True)[0][unknown_voltage]
    for idx, c in enumerate(voltages_):
        solved = solved.subs(voltage[idx + 1], c)
    return expr_to_quantity(solved)
