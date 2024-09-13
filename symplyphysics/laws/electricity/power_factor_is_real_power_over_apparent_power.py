"""
Power factor is real power over apparent power
==============================================

Power factor of an AC circuit is a dimensionless physical quantity expressing the circuit's
efficiency. It is the ratio of real power to apparent power of the circuit.

..
    TODO Move law to circuits folder?
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input,
    validate_output, symbols, clone_symbol)

apparent_power = clone_symbol(symbols.power, display_symbol="S")
"""
Apparent power of the circuit, which is the absolute value of the complex
power.
"""

real_power = clone_symbol(symbols.power, symbols="P")
"""
Real power, or active power, of the circuit.
"""

power_factor = symbols.power_factor
"""
Power factor of the circuit.
"""

law = Eq(power_factor, real_power / apparent_power)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(active_power_=real_power, full_power_=apparent_power)
@validate_output(power_factor)
def calculate_power_factor(active_power_: Quantity, full_power_: Quantity) -> Quantity:
    result_factor_expr = solve(law, power_factor, dict=True)[0][power_factor]
    result_expr = result_factor_expr.subs({real_power: active_power_, apparent_power: full_power_})
    return Quantity(result_expr)
