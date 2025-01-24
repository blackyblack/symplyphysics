"""
Mechanical energy is constant
=============================

Mechanical energy is constant in a system that has only gravitational forces or in an
otherwise idealized system — that is, one lacking dissipative forces, such as friction
and air resistance, or one in which such forces can be reasonably neglected.

**Notes:**

#. SymPy does not have a proper way to represent constant energy. We use it's derivative
   over time instead. Derivative of the constant value is zero.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Conservation_of_energy>`__.
"""

from sympy import (Eq, dsolve, Derivative)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_function)

time = symbols.time
"""
:symbols:`time`.
"""

mechanical_energy = clone_as_function(symbols.mechanical_energy, [time])
"""
:symbols:`mechanical_energy`.
"""

law = Eq(Derivative(mechanical_energy(time), time), 0)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mechanical_energy_before_=mechanical_energy)
@validate_output(mechanical_energy)
def calculate_energy_after(mechanical_energy_before_: Quantity) -> Quantity:
    solved = dsolve(law, mechanical_energy(time))
    result_expr = solved.subs("C1", mechanical_energy_before_).rhs
    return Quantity(result_expr)
