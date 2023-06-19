from typing import List
from sympy import (Eq, solve, symbols, Idx, IndexedBase, Sum, Add, Function as SymFunction)
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


def eval_sum(*args):
    n = len(args)
    return n
    #return prod(
    #    prod(args[j] - args[i] for j in range(i + 1, n))
    #    / factorial(i) for i in range(n))
    # converting factorial(i) to int is slightly faster


class MySumm(SymFunction):
    is_integer = True

    @classmethod
    def eval(cls, *args):
        return eval_sum(*args)

    def doit(self, **hints):
        return eval_sum(*self.args)


voltages = Symbol("voltages", units.voltage)
law2 = Eq(MySumm(voltages), 0)


def print() -> str:
    return print_expression(law2)


@validate_input(voltages_=units.voltage)
@validate_output(units.voltage)
def calculate_voltage(voltages_: List[Quantity]) -> Quantity:
    voltage_symbols = [Symbol("voltage_" + str(i), units.voltage) for i in range(len(voltages_) + 1)]
    unknown_voltage = voltage_symbols[len(voltages_)]
    voltages_law = law2.subs(voltages, tuple(voltage_symbols))
    solved = solve(voltages_law, unknown_voltage, dict=True)[0][unknown_voltage]
    for idx, c in enumerate(voltages_):
        solved = solved.subs(voltage[idx + 1], c)
    return expr_to_quantity(solved)
