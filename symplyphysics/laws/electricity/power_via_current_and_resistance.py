"""
Power via current and resistance
================================

Electric power can be expressed using current flowing through a
conductor and its resistance. Applied to AC circuits, the power in this law is the
instantaneous power in the circuit.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import (
    current_is_voltage_over_resistance as ohm_law,
    power_via_voltage_and_resistance as power_law,
)

power = Symbol("power", units.power)
"""
Electric power.

Symbol:
    :code:`P`
"""

current = Symbol("current", units.current)
"""
Electric current.

Symbol:
    :code:`I`
"""

resistance = Symbol("resistance", units.impedance)
r"""
Electrical resistance.

Symbol:
    :code:`R`
"""

law = Eq(power, current**2 * resistance)
r"""
:code:`P = I^2 * R`

Latex:
    .. math::
        P = I^2 R
"""

# Derive law using Ohm's law and power law

_ohm_eqn = ohm_law.law.subs({
    ohm_law.current: current,
    ohm_law.resistance: resistance,
    ohm_law.voltage: power_law.voltage,
})

_power_eqn = power_law.law.subs({
    power_law.power: power,
    power_law.resistance: resistance,
})

_power_derived = solve(
    (_ohm_eqn, _power_eqn),
    (power_law.voltage, power),
    dict=True,
)[0][power]

# Check if derived power is same as declared
assert expr_equals(_power_derived, law.rhs)


@validate_input(current_=current, resistance_=resistance)
@validate_output(power)
def calculate_heat_power(current_: Quantity, resistance_: Quantity) -> Quantity:
    result_power_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_power_expr.subs({current: current_, resistance: resistance_})
    return Quantity(result_expr)
