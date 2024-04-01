from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, SymbolIndexed, SumIndexed, global_index)

# Description
## If inductors are connected in series, total inductance is a sum of inductances of each inductor.
## Law: L_serial = sum(L[i]), where
## L_serial is total inductance,
## L[i] is inductance of i-th inductor.

# Conditions
## All inductors are NOT magnetically coupled.

serial_inductance = Symbol("serial_inductance", units.inductance)
inductance = SymbolIndexed("inductance", units.inductance)
law = Eq(serial_inductance, SumIndexed(inductance[global_index], global_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(inductances_=inductance)
@validate_output(serial_inductance)
def calculate_serial_inductance(inductances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(inductances_)))
    inductances_law = law.subs(global_index, local_index)
    inductances_law = inductances_law.doit()
    solved = solve(inductances_law, serial_inductance, dict=True)[0][serial_inductance]
    for i, v in enumerate(inductances_):
        solved = solved.subs(inductance[i + 1], v)
    return Quantity(solved)
