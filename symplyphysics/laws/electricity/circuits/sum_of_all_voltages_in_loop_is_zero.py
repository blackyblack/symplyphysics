from typing import List
from sympy import (Add, Eq, solve)
from symplyphysics import (Symbol, units, expr_to_quantity, Quantity, print_expression,
    validate_input, validate_output)

# Description
## sum(U) = 0
## Where U is a voltage on the element of an electrical loop.
## Loop is a closed sub-circuit in an electrical circuit.
## In other words, sum of all voltage sources in the loop equals to sum of all voltage consumers in this loop.
## This property of electrical loop is called Kirchhoff law #2.
## This law also demonstrates that work to move the unit of charge in electrical field along the closed path is zero.

voltages = Symbol("voltages", units.voltage)

law = Eq(Add(voltages), 0)


def print() -> str:
    return print_expression(law)


@validate_input(voltages_=units.voltage)
@validate_output(units.voltage)
def calculate_voltage(voltages_: List[Quantity]) -> Quantity:
    voltage_symbols = [
        Symbol("voltage_" + str(i), units.voltage) for i in range(0,
        len(voltages_) + 1)
    ]
    voltages_law = law.subs(voltages, tuple(voltage_symbols))
    unknown_voltage = voltage_symbols[len(voltages_)]
    solved = solve(voltages_law, unknown_voltage, dict=True)[0][unknown_voltage]
    for idx, v in enumerate(voltages_):
        solved = solved.subs(voltage_symbols[idx], v)
    return expr_to_quantity(solved)
