"""
Dynamic viscosity of gas from temperature
=========================================

Viscosity is a phenomenon of molecular transport, the property of fluid bodies (liquids and gases)
to resist the movement of one part of them relative to another. As a result, the macroscopic work
expended on this movement is dissipated as heat. Unlike liquids, the viscosity of gases increases
with increasing temperature, whereas for liquids it decreases with increasing temperature.

**Conditions:**

#. The gas is ideal.

**Links:**

#. `Sutherland model <https://en.wikipedia.org/wiki/Temperature_dependence_of_viscosity#Sutherland_model>`__.

..
    TODO Rename file
"""

from sympy import Eq, solve, Rational
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

dynamic_viscosity = symbols.dynamic_viscosity
"""
:symbols:`dynamic_viscosity` of the gas at the given :attr:`~temperature`.
"""

reference_dynamic_viscosity = clone_as_symbol(symbols.dynamic_viscosity, subscript="0")
"""
:symbols:`dynamic_viscosity` of the gas at the :attr:`~reference_temperature`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` at which the viscosity value is calculated.
"""

reference_temperature = clone_as_symbol(symbols.temperature, subscript="0")
"""
:symbols:`temperature` at which the reference viscosity value is calculated.
"""

sutherland_constant = symbols.sutherland_constant
"""
:symbols:`sutherland_constant` of the gas.
"""

law = Eq(
    dynamic_viscosity,
    reference_dynamic_viscosity * ((reference_temperature + sutherland_constant) /
    (temperature + sutherland_constant)) * (temperature / reference_temperature)**Rational(3, 2))
"""
:laws:symbol::

:laws:latex::
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
