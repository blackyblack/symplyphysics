from sympy import Eq
from symplyphysics import (
    Quantity,
    Symbol,
    units,
    validate_input,
    validate_output,
)

# Description
## The spacetime interval is invariant under Lorentz transformations, i.e. it remains the same in all
## inertial frames of reference.

# Law: (Δs)_2 = (Δs)_1
## Δs - spacetime interval
## Subscripts 2, 1 designate two inertial reference frames connected by a Lorentz transformation

# Links: Wikipedia, follows from text <https://en.wikipedia.org/wiki/Spacetime#Spacetime_interval>

first_spacetime_interval = Symbol("first_spacetime_interval", units.length)
second_spacetime_interval = Symbol("second_spacetime_interval", units.length)

law = Eq(second_spacetime_interval, first_spacetime_interval)

# TODO: derive from the definition of spacetime interval and Lorentz transformations


@validate_input(first_spacetime_interval_=first_spacetime_interval)
@validate_output(second_spacetime_interval)
def calculate_second_spacetime_interval(first_spacetime_interval_: Quantity) -> Quantity:
    result = law.rhs.subs({first_spacetime_interval: first_spacetime_interval_})
    return Quantity(result)
