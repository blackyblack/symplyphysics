from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input, validate_output)

# Description
## Resistance of the wire is proportional to it's length and resistivity and inversely proportional to it's cross-sectional area.
## Law: R = ro * l / S, where
## R is resistance of the wire
## ro is resistivity of wire material
## l is length of this wire
## S is wire cross-section.

resistance = Symbol("resistance", units.impedance)
resistivity = Symbol("resistivity", units.impedance * units.length)
wire_length = Symbol("wire_length", units.length)
cross_section = Symbol("cross_section", units.length**2)

law = Eq(resistance, resistivity * wire_length / cross_section)


#def print() -> str:
#    return print_expression(law)


@validate_input(resistivity_=resistivity, wire_length_=wire_length, cross_section_=cross_section)
#@validate_output(resistance)
def calculate_resistance(resistivity_: Quantity, wire_length_: Quantity, cross_section_: Quantity) -> Quantity:
    result_resistance_expr = solve(law, resistance, dict=True)[0][resistance]
    print(f"{result_resistance_expr}")
    result_expr = result_resistance_expr.subs({resistivity: resistivity_, wire_length: wire_length_, cross_section: cross_section_})
    print(f"{result_expr}")
    return expr_to_quantity(result_expr)


resistivity1 = Quantity(0.0172 * units.ohm * (units.milli * units.meter)**2 / units.meter)
wire_length1 = Quantity(3 * units.meter)
cross_section1 = Quantity(2 * (units.milli * units.meter)**2)

resistance = calculate_resistance(resistivity1, wire_length1, cross_section1)
print(f"{resistance}")