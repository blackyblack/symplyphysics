r"""
Critical molar volume
=====================

Critical molar volume of a van der Waals fluid is proportional to the excluded volume
parameter :math:`b` of the van der Waals equation. See :ref:`vdw_critical_parameters_def`.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

critical_molar_volume = Symbol("critical_molar_volume", units.volume / units.amount_of_substance)
r"""
Critical volume of the van der Waals fluid per unit amount of substance.

Symbol:
    :code:`v_c`

Latex:
    :math:`v_\text{c}`
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

law = Eq(critical_molar_volume, 3 * excluded_volume_parameter)
r"""
:code:`v_c = 3 * b`

Latex:
    .. math::
        v_c = 3 b
"""


@validate_input(molecules_volume_parameter_=excluded_volume_parameter)
@validate_output(critical_molar_volume)
def calculate_critical_molar_volume(molecules_volume_parameter_: Quantity) -> Quantity:
    result = law.rhs.subs(excluded_volume_parameter, molecules_volume_parameter_)
    return Quantity(result)
