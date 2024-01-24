from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohm_law

# Description
## The voltage in the circuit depends on the current in the circuit, as well as on the internal
## resistance of the circuit (resistance of current sources) and the external resistance of the circuit.
## Law is: U = I * (r + R), where
## U - voltage,
## I - current,
## r - internal resistance,
## R - external resistance.

voltage = Symbol("voltage", units.voltage)

current = Symbol("current", units.current)
inner_resistance = Symbol("inner_resistance", units.impedance)
outer_resistance = Symbol("outer_resistance", units.impedance)

law = Eq(voltage, current * (inner_resistance + outer_resistance))

# This law might be derived via Ohm's law.

ohm_law_applied = ohm_law.law.subs({
    ohm_law.voltage: ohm_law.voltage,
    ohm_law.current: current,
    ohm_law.resistance: inner_resistance + outer_resistance
})
voltage_derived = solve(ohm_law_applied, ohm_law.voltage, dict=True)[0][ohm_law.voltage]

# Check if derived voltage is same as declared.
assert expr_equals(voltage_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(current_=current, inner_resistance_=inner_resistance, outer_resistance_=outer_resistance)
@validate_output(voltage)
def calculate_voltage(current_: Quantity, inner_resistance_: Quantity, outer_resistance_: Quantity) -> Quantity:
    result_voltage_expr = solve(law, voltage, dict=True)[0][voltage]
    result_expr = result_voltage_expr.subs({current: current_, inner_resistance: inner_resistance_, outer_resistance: outer_resistance_})
    return Quantity(result_expr)
