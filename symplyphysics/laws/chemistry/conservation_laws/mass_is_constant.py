"""
Mass is constant
================

The mass is constant in a system that is closed, and mass is not transformed to energy.

**Conditions:**

#. System in a closed impenetrable volume, that is, molecules/atoms cannot leave it and
   they are always inside.
#. Mass is not transformed to energy, for example due to annihilation.

**Note:**

#. SymPy does not have a proper way to represent constant mass. We use it's derivative
   over time instead. Derivative of the constant value is zero.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Conservation_of_mass>`__.
"""

from sympy import (Eq, dsolve, Derivative)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_function)

time = symbols.time
"""
:symbols:`time`.
"""

mass = clone_as_function(symbols.mass, [time])
"""
:symbols:`mass` as a function of :attr:`~time`.
"""

law = Eq(Derivative(mass(time), time), 0)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_before_=mass)
@validate_output(mass)
def calculate_mass_after(mass_before_: Quantity) -> Quantity:
    solved = dsolve(law, mass(time))
    result_expr = solved.subs("C1", mass_before_).rhs
    return Quantity(result_expr)
