from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, print_expression, Symbol, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols
from typing import Sequence

# Description
## When the capacitors are connected in series, the voltage is equal to the sum of the voltages on the individual capacitors.
## A series connection is used to increase the breakdown voltage of the capacitors.

## Law: V = sum(V_i)
## Where:
## V is total voltage
## V_i is the voltage of one of the connected capacitors

total_voltage = Symbol("total_voltage", units.voltage)
voltages = Symbol("voltages", units.voltage)

law = Eq(total_voltage, SumArray(voltages), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(total_voltage_=total_voltage)
@validate_output(units.voltage)
def calculate_total_voltage(total_voltage_: Sequence[Quantity]) -> Quantity:
    voltages_symbols = tuple_of_symbols("voltage", units.voltage,
        len(total_voltage_))
    voltages_law = law.subs(voltages, voltages_symbols).doit()
    solved = solve(voltages_law, total_voltage, dict=True)[0][total_voltage]
    for (from_, to_) in zip(voltages_symbols, total_voltage_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
