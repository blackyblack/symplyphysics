r"""
Charge is quantized
===================

Electric charge is quantized, i.e. restricted to certain values.

**Notation:**

#. :quantity_notation:`elementary_charge`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Elementary_charge>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

charge = symbols.charge
"""
Total :symbols:`charge` of the body.
"""

integer_factor = symbols.whole_number
"""
:symbols:`whole_number`.
"""

law = Eq(charge, integer_factor * quantities.elementary_charge)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(integer_factor_=integer_factor)
@validate_output(charge)
def calculate_charge(integer_factor_: int) -> Quantity:
    charge_expr = law.rhs
    charge_value = charge_expr.subs(integer_factor, integer_factor_)
    return Quantity(charge_value)
