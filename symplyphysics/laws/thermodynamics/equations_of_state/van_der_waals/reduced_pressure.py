r"""
Reduced pressure
================

See :ref:`vdw_reduced_units_def`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#Critical_point_and_corresponding_states>`__.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    convert_to_float,
)

reduced_pressure = Symbol("reduced_pressure", dimensionless)
r"""
Reduced pressure of the van der Waals fluid.

Symbol:
    :code:`p*`

Latex:
    :math:`p^*`
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure of the van der Waals fluid.

Symbol:
    :code:`p`
"""

critical_pressure = Symbol("critical_pressure", units.pressure)
r"""
See :doc:`laws.thermodynamics.equations_of_state.van_der_waals.critical_pressure`

Symbol:
    :code:`p_c`

Latex:
    :math:`p_\text{c}`
"""

law = Eq(reduced_pressure, pressure / critical_pressure)
r"""
:code:`p* = p / p_c`

Latex:
    .. math::
        p^* = \frac{p}{p_\text{c}}
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
