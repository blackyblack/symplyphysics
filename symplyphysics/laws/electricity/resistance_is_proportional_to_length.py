from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)

# Description
## Resistance of the wire is proportional to it's length and resistivity and inversely proportional to it's section.
## Law: R = ro * l / S, where
## R is resistance of the wire
## ro is resistivity of wire material
## l is length of this wire
## S is wire section.

resistance = Symbol("resistance", units.impedance)
resistivity = Symbol("resistivity", units.impedance * units.length)
wire_length = Symbol("length", units.length)
section = Symbol("section", units.length * units.length)

law = Eq(resistance, resistivity * wire_length / section)

def print() -> str:
    return print_expression(law)

@validate_input_symbols(resistivity_=resistivity, wire_length_=wire_length, section_=section)
@validate_output_symbol(resistance)
def calculate_resistance(resistivity_: Quantity, wire_length_: Quantity, section_: Quantity) -> Quantity:
    result_resistance_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_resistance_expr.subs({resistivity: resistivity_, wire_length: wire_length_, section_: section})
    return expr_to_quantity(result_expr)
