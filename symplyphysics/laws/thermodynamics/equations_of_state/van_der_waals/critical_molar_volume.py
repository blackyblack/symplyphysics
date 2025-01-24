r"""
Critical molar volume
=====================

Critical molar volume of a van der Waals fluid is proportional to the excluded volume
parameter :math:`b` of the van der Waals equation. See :ref:`vdw_critical_parameters_def`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#Critical_point_and_corresponding_states>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

critical_molar_volume = clone_as_symbol(symbols.molar_volume, display_symbol="v_cm", display_latex="v_{\\text{c},\\text{m}}")
"""
Critical :symbols:`molar_volume` of the van der Waals fluid.
"""

excluded_volume_parameter = symbols.excluded_volume_parameter
"""
:symbols:`excluded_volume_parameter`.
"""

law = Eq(critical_molar_volume, 3 * excluded_volume_parameter)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(molecules_volume_parameter_=excluded_volume_parameter)
@validate_output(critical_molar_volume)
def calculate_critical_molar_volume(molecules_volume_parameter_: Quantity) -> Quantity:
    result = law.rhs.subs(excluded_volume_parameter, molecules_volume_parameter_)
    return Quantity(result)
