r"""
Charge is quantized
===================

Electric charge is quantized, i.e. restricted to certain values.

**Notation:**

#. :math:`e` is the elementary charge.
"""

from sympy import Eq
from symplyphysics import (
    SymbolNew,
    dimensionless,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

charge = symbols.charge
"""
Total charge of the body.
"""

integer_factor = SymbolNew("n", dimensionless)  # positive, zero or negative
"""
Integer factor of any sign.
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
