from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## If resistors are connected in series, total resistance is a sum of resistances of each resistor.
## Law: R_serial = sum(R[i]), where
## R_serial is total resistance,
## R[i] is resistance of i-th resistor.

# Conditions:
# - it is valid for non-alternating current, or when there is no reactive load.

resistances = Symbol("resistances", units.impedance)
serial_resistance = Symbol("serial_resistance", units.impedance)
law = Eq(serial_resistance, SumArray(resistances), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(resistances_=resistances)
@validate_output(units.impedance)
def calculate_serial_resistance(resistances_: list[Quantity]) -> Quantity:
    resistance_symbols = tuple_of_symbols("resistance", units.impedance, len(resistances_))
    resistances_law = law.subs(resistances, resistance_symbols).doit()
    solved = solve(resistances_law, serial_resistance, dict=True)[0][serial_resistance]
    for (from_, to_) in zip(resistance_symbols, resistances_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
