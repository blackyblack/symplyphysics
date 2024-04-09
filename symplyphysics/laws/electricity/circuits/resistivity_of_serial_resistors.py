from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, SymbolIndexed, SumIndexed, global_index)

# Description
## If resistors are connected in series, total resistance is a sum of resistances of each resistor.
## Law: R_serial = sum(R[i]), where
## R_serial is total resistance,
## R[i] is resistance of i-th resistor.

# Conditions:
# - it is valid for non-alternating current, or when there is no reactive load.

resistances = Symbol("resistances", units.impedance)
serial_resistance = Symbol("serial_resistance", units.impedance)
resistance = SymbolIndexed("resistance", units.impedance)
law = Eq(serial_resistance, SumIndexed(resistance[global_index], global_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(resistances_=resistance)
@validate_output(serial_resistance)
def calculate_serial_resistance(resistances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(resistances_)))
    resistances_law = law.subs(global_index, local_index)
    resistances_law = resistances_law.doit()
    solved = solve(resistances_law, serial_resistance, dict=True)[0][serial_resistance]
    for i, v in enumerate(resistances_):
        solved = solved.subs(resistance[i + 1], v)
    return Quantity(solved)
