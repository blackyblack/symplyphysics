from typing import List
from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, print_expression,
    Symbol, validate_input, validate_output)
from symplyphysics.core.operations.sum_array import SumArray

# Description
## sum(U) = 0
## Where U is a voltage on the element of an electrical loop.
## Loop is a closed sub-circuit in an electrical circuit.
## In other words, sum of all voltage sources in the loop equals to sum of all voltage consumers in this loop.
## This property of electrical loop is called Kirchhoff law #2.
## This law also demonstrates that work to move the unit of charge in electrical field along the closed path is zero.

voltages = Symbol("voltages", units.voltage)
law = Eq(SumArray(voltages), 0, evaluate=False)


def print() -> str:
    return print_expression(law)


@validate_input(voltages_=voltages)
@validate_output(units.voltage)
def calculate_voltage(voltages_: List[Quantity]) -> Quantity:
    voltage_symbols = tuple(Symbol("voltage_" + str(i), units.voltage) for i in range(len(voltages_) + 1))
    unknown_voltage = voltage_symbols[len(voltages_)]
    voltages_law = law.subs(voltages, tuple(voltage_symbols)).doit()
    solved = solve(voltages_law, unknown_voltage, dict=True)[0][unknown_voltage]
    for idx in range(len(voltage_symbols) - 1):
        solved = solved.subs(voltage_symbols[idx], voltages_[idx])
    return expr_to_quantity(solved)
