"""
Spacetime interval is Lorentz invariant
=======================================

The spacetime interval is invariant under Lorentz transformations, i.e. it remains the same in all
inertial frames of reference.

**Links:**

#. `Wikipedia, follows from text <https://en.wikipedia.org/wiki/Spacetime#Spacetime_interval>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

first_spacetime_interval = clone_as_symbol(symbols.spacetime_interval, subscript="1")
"""
:symbols:`spacetime_interval` in one inertial frame of reference.
"""

second_spacetime_interval = clone_as_symbol(symbols.spacetime_interval, subscript="2")
"""
:symbols:`spacetime_interval` in another inertial frame of reference related to the former
by a Lorentz transformation.
"""

law = Eq(second_spacetime_interval, first_spacetime_interval)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive from the definition of spacetime interval and Lorentz transformations


@validate_input(first_spacetime_interval_=first_spacetime_interval)
@validate_output(second_spacetime_interval)
def calculate_second_spacetime_interval(first_spacetime_interval_: Quantity) -> Quantity:
    result = law.rhs.subs({first_spacetime_interval: first_spacetime_interval_})
    return Quantity(result)
