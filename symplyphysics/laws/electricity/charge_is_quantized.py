r"""
Charge is quantized
===================

Electric charge is quantized, i.e. restricted to certain values.

**Notation:**

#. :math:`e` is the elementary charge.
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

charge = Symbol("charge", units.charge)
"""
Total charge of the body.

Symbol:
    :code:`q`
"""

integer_factor = Symbol("factor", dimensionless)  # positive, zero or negative
"""
Integer factor of any sign.

Symbol:
    :code:`n`
"""

law = Eq(charge, integer_factor * units.elementary_charge)
r"""
:code:`q = n * e`

Latex:
    .. math::
        q = n e
"""


@validate_input(integer_factor_=integer_factor)
@validate_output(charge)
def calculate_charge(integer_factor_: int) -> Quantity:
    charge_expr = law.rhs
    charge_value = charge_expr.subs(integer_factor, integer_factor_)
    return Quantity(charge_value)
