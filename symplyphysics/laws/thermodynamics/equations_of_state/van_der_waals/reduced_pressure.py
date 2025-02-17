r"""
Reduced pressure
================

See :ref:`vdw_reduced_units_def`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#Critical_point_and_corresponding_states>`__.
"""

from sympy import Eq
from symplyphysics import (
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

reduced_pressure = Symbol("p_r", dimensionless)
"""
Reduced :symbols:`pressure` of the van der Waals fluid.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` of the van der Waals fluid.
"""

critical_pressure = clone_as_symbol(symbols.pressure, display_symbol="p_c", display_latex="p_\\text{c}")
r"""
Critical :symbols:`pressure`. See :doc:`laws.thermodynamics.equations_of_state.van_der_waals.critical_pressure`.
"""

law = Eq(reduced_pressure, pressure / critical_pressure)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    pressure_=pressure,
    critical_pressure_=critical_pressure,
)
@validate_output(reduced_pressure)
def calculate_reduced_pressure(
    pressure_: Quantity,
    critical_pressure_: Quantity,
) -> float:
    result = law.rhs.subs({
        pressure: pressure_,
        critical_pressure: critical_pressure_,
    })
    return convert_to_float(result)
