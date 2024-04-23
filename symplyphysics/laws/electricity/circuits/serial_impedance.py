from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, SymbolIndexed, SumIndexed, global_index)

# Description
## If elements of circuit are connected in series, total impedance is a sum of impedances of each element.
## Law: Z_serial = sum(Z[i]), where
## Z_serial is total impedance,
## Z[i] is impedance of i-th element of circuit.

serial_impedance = Symbol("serial_impedance", units.impedance)
impedance = SymbolIndexed("impedance", units.impedance)
law = Eq(serial_impedance, SumIndexed(impedance[global_index], global_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(impedances_=impedance)
@validate_output(serial_impedance)
def calculate_serial_impedance(impedances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(impedances_)))
    impedances_law = law.subs(global_index, local_index)
    impedances_law = impedances_law.doit()
    solved = solve(impedances_law, serial_impedance, dict=True)[0][serial_impedance]
    for i, v in enumerate(impedances_):
        solved = solved.subs(impedance[i + 1], v)
    return Quantity(solved)
