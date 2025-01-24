r"""
Critical temperature
====================

Critical temperature of a van der Waals fluid depends on the parameters :math:`a` and :math:`b`
of the van der Waals equation and the molar gas constant :math:`R`. See :ref:`vdw_critical_parameters_def`.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#Critical_point_and_corresponding_states>`__.
"""

from sympy import Eq
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
)

critical_temperature = clone_as_symbol(symbols.temperature, subscript="\\text{c}")
"""
Critical :symbols:`temperature` of the van der Waals fluid.
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
    critical_temperature,
    (8 * attractive_forces_parameter) /
    (27 * quantities.molar_gas_constant * excluded_volume_parameter),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    bonding_forces_parameter_=attractive_forces_parameter,
    molecules_volume_parameter_=excluded_volume_parameter,
)
@validate_output(critical_temperature)
def calculate_critical_temperature(
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
) -> Quantity:
    temperature_ = law.rhs.subs({
        attractive_forces_parameter: bonding_forces_parameter_,
        excluded_volume_parameter: molecules_volume_parameter_,
    })
    return Quantity(temperature_)
