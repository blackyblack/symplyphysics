"""
Reduced temperature
===================

See :ref:`vdw_reduced_units_def`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#Critical_point_and_corresponding_states>`__.
"""

from sympy import Eq
from symplyphysics import (
    clone_as_symbol,
    symbols,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    convert_to_float,
)

reduced_temperature = Symbol("T_r", dimensionless)
"""
Reduced :symbols:`temperature` of the van der Waals fluid.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the van der Waals fluid.
"""

critical_temperature = clone_as_symbol(symbols.temperature, display_symbol="T_c", display_latex="T_\\text{c}")
"""
See :doc:`laws.thermodynamics.equations_of_state.van_der_waals.critical_temperature`
and :symbols:`temperature`.
"""

law = Eq(reduced_temperature, temperature / critical_temperature)
"""
:laws:symbol::

:laws:latex::
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
