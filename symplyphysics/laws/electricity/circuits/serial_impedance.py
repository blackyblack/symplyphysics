from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.operations.sum_array import SumArray
from symplyphysics.core.symbols.symbols import tuple_of_symbols

# Description
## If elements of circuit are connected in series, total impedance is a sum of impedances of each element.
## Law: Z_serial = sum(Z[i]), where
## Z_serial is total impedance,
## Z[i] is impedance of i-th element of circuit.

impedances = Symbol("impedance", units.impedance)
serial_impedance = Symbol("serial_impedance", units.impedance)
law = Eq(serial_impedance, SumArray(impedances), evaluate=False)


def print_law() -> str:
    return print_expression(law)


@validate_input(impedance_=impedances)
@validate_output(units.impedance)
def calculate_serial_impedance(impedance_: list[Quantity]) -> Quantity:
    impedance_symbols = tuple_of_symbols("impedance", units.impedance, len(impedance_))
    impedance_law = law.subs(impedances, impedance_symbols).doit()
    solved = solve(impedance_law, serial_impedance, dict=True)[0][serial_impedance]
    for (from_, to_) in zip(impedance_symbols, impedance_):
        solved = solved.subs(from_, to_)
    return Quantity(solved)
