from typing import Sequence
from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, print_expression, validate_input, validate_output,
    SymbolIndexed, SumIndexed, global_index)

# Description
## sum(U) = 0
## Where U is a voltage on the element of an electrical loop.
## Loop is a closed sub-circuit in an electrical circuit.
## In other words, sum of all voltage sources in the loop equals to sum of all voltage consumers in this loop.
## This property of electrical loop is also known as second Kirchhoff law.
## This law also demonstrates that work to move the unit of charge in electrical field along the closed path is zero.

voltage = SymbolIndexed("voltage", units.voltage)
law = Eq(SumIndexed(voltage[global_index], global_index), 0)


def print_law() -> str:
    return print_expression(law)


@validate_input(voltages_=voltage)
@validate_output(units.voltage)
def calculate_voltage(voltages_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(voltages_) + 1))
    voltages_law = law.subs(global_index, local_index)
    voltages_law = voltages_law.doit()
    unknown_voltage = voltage[len(voltages_) + 1]
    solved = solve(voltages_law, unknown_voltage, dict=True)[0][unknown_voltage]
    for i, v in enumerate(voltages_):
        solved = solved.subs(voltage[i + 1], v)
    return Quantity(solved)
