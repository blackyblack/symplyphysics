from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## If inductors are connected in series, total inductance is a sum of inductances of each inductor.
## Law: L_serial = sum(L[i]), where
## L_serial is total inductance,
## L[i] is inductance of i-th inductor.

# Conditions
## All inductors are NOT magnetically coupled.

inductances = Symbol("inductances", units.inductance)
serial_inductance = Symbol("serial_inductance", units.inductance)
law = Eq(serial_inductance, SumArray(inductances), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(inductances_=inductances)
@validate_output(units.inductance)
def calculate_serial_inductance(inductances_: list[Quantity]) -> Quantity:
    inductance_symbols = tuple_of_symbols("inductance", units.inductance, len(inductances_))
    inductances_law = law.subs(inductances, inductance_symbols).doit()
    solved = solve(inductances_law, serial_inductance, dict=True)[0][serial_inductance]
    for (from_, to_) in zip(inductance_symbols, inductances_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
