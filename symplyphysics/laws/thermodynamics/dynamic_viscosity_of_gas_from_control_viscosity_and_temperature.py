"""
Dynamic viscosity of gas from temperature
=========================================

Viscosity is a phenomenon of molecular transport, the property of fluid bodies (liquids and gases)
to resist the movement of one part of them relative to another. As a result, the macroscopic work
expended on this movement is dissipated as heat. Unlike liquids, the viscosity of gases increases
with increasing temperature, whereas for liquids it decreases with increasing temperature.

**Conditions:**

#. The gas is ideal.

..
    TODO Rename file
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_as_symbol, symbols, units, Quantity, Symbol, validate_input,
    validate_output)

dynamic_viscosity = Symbol("dynamic_viscosity", units.pressure * units.time)
r"""
Dynamic viscosity of the gas at temperature :math:`T`.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

reference_dynamic_viscosity = Symbol("reference_dynamic_viscosity", units.pressure * units.time)
r"""
Dynamic viscosity of the gas at the reference temperature :math:`T_0`.

Symbol:
    :code:`mu_0`

Latex:
    :math:`\mu_0`
"""

temperature = symbols.temperature
"""
:symbols:`temperature` at which the viscosity value is calculated.
"""

reference_temperature = clone_as_symbol(symbols.temperature,
    display_symbol="T0",
    display_latex="T_0")
"""
:symbols:`temperature` at which the reference viscosity value is calculated.
"""

sutherland_constant = Symbol("sutherland_constant", units.temperature)
"""
Sutherland constant of the gas.

Symbol:
    :code:`C`
"""

law = Eq(
    dynamic_viscosity,
    reference_dynamic_viscosity * ((reference_temperature + sutherland_constant) /
    (temperature + sutherland_constant)) * (temperature / reference_temperature)**1.5)
r"""
:code:`mu = mu_0 * ((T0 + C) / (T + C)) * (T / T0)^(3/2)`

Latex:
    .. math::
        \mu = \mu_0 \left( \frac{T_0 + C}{T + C} \right) \left( \frac{T}{T_0} \right)^{3 / 2}
"""


@validate_input(control_viscosity_=reference_dynamic_viscosity,
    control_temperature_=reference_temperature,
    sutherland_constant_=sutherland_constant,
    set_temperature_=temperature)
@validate_output(dynamic_viscosity)
def calculate_dynamic_viscosity(control_viscosity_: Quantity, control_temperature_: Quantity,
    sutherland_constant_: Quantity, set_temperature_: Quantity) -> Quantity:
    result_expr = solve(law, dynamic_viscosity, dict=True)[0][dynamic_viscosity]
    result_dynamic_viscosity = result_expr.subs({
        reference_dynamic_viscosity: control_viscosity_,
        reference_temperature: control_temperature_,
        sutherland_constant: sutherland_constant_,
        temperature: set_temperature_
    })
    return Quantity(result_dynamic_viscosity)
