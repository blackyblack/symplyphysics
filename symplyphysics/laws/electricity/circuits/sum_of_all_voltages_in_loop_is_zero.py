from functools import reduce
import operator
from typing import Iterable, List
from sympy import (Eq, solve, symbols, Idx, IndexedBase, Sum, Expr, sympify)
from symplyphysics import (units, expr_to_quantity, Quantity, print_expression,
    Symbol, validate_input, validate_output)

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


class SumArray(Expr):
    """
    Represents unevaluated Sum over array.

    """

    def __new__(cls, *array):
        obj = Expr.__new__(cls, *array)
        return obj

    def doit(self, **hints):
        return reduce(operator.add, self._args) if isinstance(self._args, Iterable) else self._args



voltages = Symbol("voltages", units.voltage)
law2 = Eq(SumArray(voltages, voltages), 0, evaluate=False)


def print() -> str:
    return print_expression(law2)


@validate_input(voltages_=units.voltage)
@validate_output(units.voltage)
def calculate_voltage(voltages_: List[Quantity]) -> Quantity:
    voltage_symbols = [Symbol("voltage_" + str(i), units.voltage) for i in range(len(voltages_) + 1)]
    sympy_voltages = sympify(voltage_symbols, strict=False)
    unknown_voltage = voltage_symbols[len(voltages_)]
    voltages_law = law2.subs(voltages, sympy_voltages)
    #voltages_law = law2
    #voltages_law = law2.subs(voltages, voltage_symbols[0])
    voltages_law = voltages_law.doit()
    solved = solve(voltages_law, unknown_voltage, dict=True)[0][unknown_voltage]
    for idx in range(len(voltage_symbols) - 1):
        solved = solved.subs(voltage_symbols[idx], voltages_[idx])
    return expr_to_quantity(solved)
