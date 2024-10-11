r"""
Critical temperature
====================

Critical temperature of a van der Waals fluid depends on the parameters :math:`a` and :math:`b`
of the van der Waals equation and the molar gas constant :math:`R`. See :ref:`vdw_critical_parameters_def`.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.
"""

from sympy import Eq
from symplyphysics import (
    clone_as_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    quantities,
)

critical_temperature = clone_as_symbol(symbols.temperature,
    display_symbol="T_c",
    display_latex="T_\\text{c}")
"""
Critical :symbols:`temperature` of the van der Waals fluid.
"""

attractive_forces_parameter = Symbol(
    "attractive_forces_parameter",
    units.pressure * (units.volume / units.amount_of_substance)**2,
)
"""
Parameter of the van der Waals equation denoting the magnitude of attractive
forces between gas molecules.

Symbol:
    :code:`a`
"""

excluded_volume_parameter = Symbol(
    "excluded_volume_parameter",
    units.volume / units.amount_of_substance,
)
"""
Parameter of the van der Waals equation denoting an excluded molar volume
due to a finite size of molecules.

Symbol:
    :code:`b`
"""

law = Eq(
    critical_temperature,
    (8 * attractive_forces_parameter) /
    (27 * quantities.molar_gas_constant * excluded_volume_parameter),
)
r"""
:code:`T_c = 8 * a / (27 * R * b)`

Latex:
    .. math::
        T_\text{c} = \frac{8 a}{27 R b}
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
