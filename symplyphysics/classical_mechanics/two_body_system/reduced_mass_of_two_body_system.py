"""
Reduced mass of a two-body system
=================================

Reduced mass is effective inertial mass in a system with two or more particles when they
are interacting with each other. This allows the two-body problem to be solved as if it
were a one-body problem.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Reduced_mass#Equation>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol

reduced_mass = clone_as_symbol(symbols.mass, display_symbol="mu", display_latex="\\mu")
"""
The reduced :symbols:`mass` of the system.
"""

first_mass = clone_as_symbol(symbols.mass, subscript="1")
"""
The :symbols:`mass` of the first body.
"""

second_mass = clone_as_symbol(symbols.mass, subscript="2")
"""
The :symbols:`mass` of the second body.
"""

law = Eq(reduced_mass, 1 / (1 / first_mass + 1 / second_mass))
"""
:laws:symbol::

:laws:latex::
"""

# This law is more of a definition of reduced mass, so it's not derivable.


@validate_input(
    first_mass_=first_mass,
    second_mass_=second_mass,
)
@validate_output(reduced_mass)
def calculate_reduced_mass(
    first_mass_: Quantity,
    second_mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        first_mass: first_mass_,
        second_mass: second_mass_,
    })
    return Quantity(result)
