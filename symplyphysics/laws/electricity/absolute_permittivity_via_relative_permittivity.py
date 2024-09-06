r"""
Absolute permittivity via relative permittivity
===============================================

Absolute permittivity of a medium can be expressed as the product of its relative permittivity
and vacuum permittivity.

**Notation:**

#. :math:`\varepsilon_0` (:code:`epsilon_0`) is **vacuum permittivity**, also called **permittivity
   of free space** or the **electric constant**.
"""

from sympy import Eq
from symplyphysics import (
    Symbol,
    units,
    dimensionless,
    Quantity,
    validate_input,
    validate_output,
)

absolute_permittivity = Symbol("absolute_permittivity", units.capacitance / units.length)
r"""
Absolute permittivity.

Symbol:
    :code:`epsilon`

Latex:
    :math:`\varepsilon`
"""

relative_permittivity = Symbol("relative_permittivity", dimensionless)
r"""
Permittivity relative to that of vacuum.

Symbol:
    :code:`epsilon_r`

Latex:
    :math:`\varepsilon_r`
"""

law = Eq(absolute_permittivity, units.vacuum_permittivity * relative_permittivity)
r"""
:code:`epsilon = epsilon_0 * epsilon_r`

Latex:
    .. math::
        \varepsilon = \varepsilon_0 \varepsilon_r
"""


@validate_input(relative_permittivity_=relative_permittivity)
@validate_output(absolute_permittivity)
def calculate_absolute_permittivity(relative_permittivity_: float) -> Quantity:
    result = law.rhs.subs(relative_permittivity, relative_permittivity_)
    return Quantity(result)
