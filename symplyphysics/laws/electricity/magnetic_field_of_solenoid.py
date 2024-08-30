"""
Magnetic field of solenoid
==========================

Near the center of the solenoid, the magnetic field is quite uniform and directly
proportional to the current in the solenoid's wire.

**Notation:**

#. :math:`\mu_0` (:code:`mu_0`) is vacuum permeability.

**Conditions:**

#. The medium is vacuum.
#. The magnetic field is measured near the center of the solenoid.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, dimensionless)

magnetic_field = Symbol("magnetic_field", units.magnetic_flux_density)
"""
Value of the magnetic field.

Symbol:
    :code:`B`
"""

current = Symbol("current", units.current)
"""
Current flowing through the solenoid.

Symbol:
    :code:`I`
"""

length = Symbol("length", units.length)
"""
Length of the solenoid.

Symbol:
    :code:`l`
"""

coil_turn_count = Symbol("coil_turn_count", dimensionless)
"""
Number of turns in the coil.

Symbol:
    :code:`N`
"""

law = Eq(magnetic_field, units.vacuum_permeability * current * coil_turn_count / length)
r"""
:code:`B = mu_0 * I * N / l`

Latex:
    .. math::
        B = \mu_0 I \frac{N}{l}
"""


@validate_input(current_=current,
    length_=length,
    number_turns_=coil_turn_count)
@validate_output(magnetic_field)
def calculate_induction(current_: Quantity, length_: Quantity,
    number_turns_: float) -> Quantity:
    if number_turns_ < 0:
        raise ValueError("Number of turns cannot be negative")
    result_expr = solve(law, magnetic_field, dict=True)[0][magnetic_field]
    result_expr = result_expr.subs({
        current: current_,
        length: length_,
        coil_turn_count: number_turns_,
    })
    return Quantity(result_expr)
