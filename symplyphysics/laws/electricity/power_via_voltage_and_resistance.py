"""
Power via voltage and resistance
================================

Power can be found using voltage and resistance.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import (
    power_via_voltage_and_current as _power_law,
    current_is_voltage_over_resistance as _ohm_law,
)

power = Symbol("power", units.power)
"""
Power.

Symbol:
    :code:`P`
"""

voltage = Symbol("voltage", units.voltage)
"""
Voltage.

Symbol:
    :code:`V`
"""

resistance = Symbol("resistance", units.impedance)
"""
Resistance.

Symbol:
    :code:`R`
"""

law = Eq(power, voltage**2 / resistance)
r"""
:code:`P = V^2 / R`

Latex:
    .. math::
        P = \frac{V^2}{R}
"""

# Derive from power and Ohm's laws

_current = _ohm_law.current

_power_eqn = _power_law.law.subs({
    _power_law.power: power,
    _power_law.current: _current,
    _power_law.voltage: voltage,
})

_ohm_eqn = _ohm_law.law.subs({
    _ohm_law.voltage: voltage,
    _ohm_law.resistance: resistance,
})

_power_expr = solve(
    (_power_eqn, _ohm_eqn),
    (power, _current),
    dict=True,
)[0][power]

assert expr_equals(_power_expr, law.rhs)


@validate_input(
    voltage_=voltage,
    resistance_=resistance,
)
@validate_output(power)
def calculate_power(
    voltage_: Quantity,
    resistance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        voltage: voltage_,
        resistance: resistance_,
    })
    return Quantity(result)
