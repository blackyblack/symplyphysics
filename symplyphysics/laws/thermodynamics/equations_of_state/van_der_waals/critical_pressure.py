r"""
Critical pressure
=================

Critical pressure in a van der Waals fluid depends on the parameters :math:`a` and
:math:`b` of the van der Waals equation. See :ref:`vdw_critical_parameters_def`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#Critical_point_and_corresponding_states>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
    symbols,
)

critical_pressure = clone_as_symbol(symbols.pressure, subscript="\\text{c}")
"""
Critical :symbols:`pressure` of the van der Waals fluid.
"""

attractive_forces_parameter = symbols.attractive_forces_parameter
"""
:symbols:`attractive_forces_parameter`.
"""

excluded_volume_parameter = symbols.excluded_volume_parameter
"""
:symbols:`excluded_volume_parameter`.
"""

law = Eq(
    critical_pressure,
    attractive_forces_parameter / (27 * excluded_volume_parameter**2),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    bonding_forces_parameter_=attractive_forces_parameter,
    molecules_volume_parameter_=excluded_volume_parameter,
)
@validate_output(critical_pressure)
def calculate_critical_pressure(
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
) -> Quantity:
    pressure_ = law.rhs.subs({
        attractive_forces_parameter: bonding_forces_parameter_,
        excluded_volume_parameter: molecules_volume_parameter_,
    })
    return Quantity(pressure_)
