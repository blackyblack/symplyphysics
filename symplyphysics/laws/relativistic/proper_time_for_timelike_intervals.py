from sympy import Eq
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    assert_equal,
)

# Description
## In relativity, proper time along a timelike world line is defined as the time as measured by a clock following
## that line. The proper time interval between two events on a world line is the change in proper time, which is 
## independent of coordinates, and is a Lorentz scalar.

# Law: Δτ = Δs / c
## Δτ - proper time
## Δs - spacetime interval between two events
## c - speed of light

# Conditions
## - The interval is timelike, i.e. `Δs` is real

proper_time = Symbol("proper_time", units.time, real=True)
spacetime_interval = Symbol("spacetime_interval", units.length)

law = Eq(proper_time, spacetime_interval / speed_of_light)


@validate_input(spacetime_interval_=spacetime_interval)
@validate_output(proper_time)
def calculate_proper_time(
    spacetime_interval_: Quantity,
) -> Quantity:
    _, imag = spacetime_interval_.scale_factor.as_real_imag()
    try:
        assert_equal(imag, 0)
    except AssertionError as e:
        raise ValueError("Interval must be real") from e

    result = law.rhs.subs({spacetime_interval: spacetime_interval_})
    return Quantity(result)
