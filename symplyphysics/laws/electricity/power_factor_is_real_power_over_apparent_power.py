"""
Power factor is real power over apparent power
==============================================

Power factor of an AC circuit is a dimensionless physical quantity expressing the circuit's
efficiency. It is the ratio of real power to apparent power of the circuit.

..
    TODO Move law to circuits folder?
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, dimensionless)

apparent_power = Symbol("apparent_power", units.power)
"""
Apparent power of the circuit, which is the absolute value of the complex
power.

Symbol:
    :code:`S`
"""

real_power = Symbol("real_power", units.power)
"""
Real power, or active power, of the circuit.

Symbol:
    :code:`P`
"""

power_factor = Symbol("power_factor", dimensionless)
r"""
Power factor of the circuit.

Symbol:
    :code:`pf`

Latex:
    :math:`\mathrm{pf}`
"""

law = Eq(power_factor, real_power / apparent_power)
r"""
:code:`pf = P / S`

Latex:
    .. math::
        \mathrm{pf} = \frac{P}{S}
"""


@validate_input(active_power_=real_power, full_power_=apparent_power)
@validate_output(power_factor)
def calculate_power_factor(active_power_: Quantity, full_power_: Quantity) -> Quantity:
    result_factor_expr = solve(law, power_factor, dict=True)[0][power_factor]
    result_expr = result_factor_expr.subs({real_power: active_power_, apparent_power: full_power_})
    return Quantity(result_expr)
