from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import (
    current_is_voltage_over_resistance as ohm_law,
    power_via_voltage_and_resistance as power_law,
)

# Description
# Dissipated heat power is proportional to current square and resistance
# P = I**2 * R
# where :
# P - is dissipated heat power of the element
# I is current flowing through this element
# R is impedance of this element

heat_power = Symbol("heat_power", units.power)
current = Symbol("current", units.current)
resistance = Symbol("resistance", units.impedance)

law = Eq(heat_power, current**2 * resistance)

# Derive law using Ohm's law and power law

_ohm_eqn = ohm_law.law.subs({
    ohm_law.current: current,
    ohm_law.resistance: resistance,
    ohm_law.voltage: power_law.voltage,
})

_power_eqn = power_law.law.subs({
    power_law.power: heat_power,
    power_law.resistance: resistance,
})

_power_derived = solve(
    (_ohm_eqn, _power_eqn),
    (power_law.voltage, heat_power),
    dict=True,
)[0][heat_power]

# Check if derived power is same as declared
assert expr_equals(_power_derived, law.rhs)


@validate_input(current_=current, resistance_=resistance)
@validate_output(heat_power)
def calculate_heat_power(current_: Quantity, resistance_: Quantity) -> Quantity:
    result_power_expr = solve(law, heat_power, dict=True)[0][heat_power]
    result_expr = result_power_expr.subs({current: current_, resistance: resistance_})
    return Quantity(result_expr)
