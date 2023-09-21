from sympy import (I, Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, angle_type)

# Description
## While the serial resistance of ideal capacitor is zero, its reactance depends on its capacitance and frequency.
## Law: Zc = -j/wС, where
## Xc is capacitor impedance,
## j is imaginary number,
## w is circular frequency,
## C is capacitance of capacitor.

capacitor_impedance = Symbol("capacitor_impedance", units.impedance)
circular_frequency = Symbol("circular_frequency", angle_type / units.time)
capacitor_capacitance = Symbol("capacitor_capacitance", units.capacitance)

law = Eq(capacitor_impedance, -I / (circular_frequency * capacitor_capacitance))

def print_law() -> str:
    return print_expression(law)


@validate_input(capacitance_=capacitor_capacitance, circular_frequency_=circular_frequency)
@validate_output(capacitor_impedance)
def calculate_impedance(capacitance_: Quantity, circular_frequency_: Quantity) -> Quantity:
    result_impedance_expr = solve(law, capacitor_impedance, dict=True)[0][capacitor_impedance]
    result_expr = result_impedance_expr.subs({
        capacitor_capacitance: capacitance_,
        circular_frequency: circular_frequency_
    })
    return Quantity(result_expr)
