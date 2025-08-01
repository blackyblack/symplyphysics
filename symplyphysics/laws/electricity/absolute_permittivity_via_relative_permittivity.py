"""
Absolute permittivity via relative permittivity
===============================================

Absolute permittivity of a medium can be expressed as the product of its relative permittivity
and vacuum permittivity.

**Notes:**

#. This equation can be seen as the definition of relative permittivity.

**Notation:**

#. :quantity_notation:`vacuum_permittivity`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Relative_permittivity#Definition>`__.
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
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity)
@validate_output(absolute_permittivity)
def calculate_absolute_permittivity(relative_permittivity_: float) -> Quantity:
    result = law.rhs.subs(relative_permittivity, relative_permittivity_)
    return Quantity(result)
