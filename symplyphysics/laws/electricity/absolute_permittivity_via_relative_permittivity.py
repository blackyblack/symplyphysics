r"""
Absolute permittivity via relative permittivity
===============================================

Absolute permittivity of a medium can be expressed as the product of its relative permittivity
and vacuum permittivity.

**Notation:**

#. :quantity_notation:`vacuum_permittivity`.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

absolute_permittivity = symbols.absolute_permittivity
"""
:symbols:`absolute_permittivity`.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity`.
"""

law = Eq(absolute_permittivity, quantities.vacuum_permittivity * relative_permittivity)
"""
:laws:symbols::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity)
@validate_output(absolute_permittivity)
def calculate_absolute_permittivity(relative_permittivity_: float) -> Quantity:
    result = law.rhs.subs(relative_permittivity, relative_permittivity_)
    return Quantity(result)
