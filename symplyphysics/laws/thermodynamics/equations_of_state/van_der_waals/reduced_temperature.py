"""
Reduced temperature
===================

Reduced units are used in the dimensionless van der Waals equation of state.
"""

from sympy import Eq
from symplyphysics import (
    clone_symbol,
    symbols,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    convert_to_float,
)

reduced_temperature = Symbol("reduced_temperature", dimensionless)
r"""
Reduced temperature of the van der Waals fluid.

Symbol:
    :code:`T*`

Latex:
    :math:`T^*`
"""

temperature = symbols.thermodynamics.temperature
r"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the van der Waals fluid.

Symbol: 
    :code:`T`
"""

critical_temperature = clone_symbol(symbols.thermodynamics.temperature, "critical_temperature")
r"""
See :doc:`laws.thermodynamics.equations_of_state.van_der_waals.critical_temperature`.

Symbol:
    :code:`T_c`

Latex:
    :math:`T_\text{c}`
"""

law = Eq(reduced_temperature, temperature / critical_temperature)
r"""
:code:`T* = T / T_c`

Latex:
    .. math::
        T^* = \frac{T}{T_\text{c}}
"""


@validate_input(
    temperature_=temperature,
    critical_temperature_=critical_temperature,
)
@validate_output(reduced_temperature)
def calculate_reduced_temperature(
    temperature_: Quantity,
    critical_temperature_: Quantity,
) -> float:
    result = law.rhs.subs({
        temperature: temperature_,
        critical_temperature: critical_temperature_,
    })
    return convert_to_float(result)
