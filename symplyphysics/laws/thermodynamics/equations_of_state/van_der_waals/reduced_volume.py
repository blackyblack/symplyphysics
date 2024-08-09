r"""
Reduced volume
==============

Reduced units are used in the dimensionless van der Waals equation of state.

**Note:**

#. Specific or molar volumes can be used in the right-hand side of the law.
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

reduced_volume = Symbol("reduced_volume", dimensionless)
r"""
Reduced volume of the van der Waals fluid.

Symbol:
    :code:`V*`

Latex:
    :math:`V^*`
"""

volume = Symbol("volume", units.volume)
"""
Volume of the van der Waals fluid.

Symbol:
    :code:`V`
"""

critical_volume = Symbol("critical_volume", units.volume)
r"""
See :doc:`laws.thermodynamics.equations_of_state.van_der_waals.critical_molar_volume`.

Symbol:
    :code:`V_c`

Latex:
    :math:`V_\text{c}`
"""

law = Eq(reduced_volume, volume / critical_volume)
r"""
:code:`V* = V / V_c`

Latex:
    .. math::
        V^* = \frac{V}{V_\text{c}}
"""


@validate_input(
    volume_=volume,
    critical_volume_=critical_volume,
)
@validate_output(reduced_volume)
def calculate_reduced_volume(
    volume_: Quantity,
    critical_volume_: Quantity,
) -> float:
    result = law.rhs.subs({
        volume: volume_,
        critical_volume: critical_volume_,
    })
    return convert_to_float(result)
