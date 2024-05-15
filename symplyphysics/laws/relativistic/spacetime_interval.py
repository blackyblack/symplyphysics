from sympy import Eq, solve
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The spacetime interval is a combination of distance and time that is invariant under Lorentz tranformations.
## It is invariant in the sense that it has the same value for all observers in any inertial reference frame.

# Law: Δs**2 = (c * Δt)**2 - Δd**2
## Δs - spacetime interval (may be complex)
## Δt - time separation of events
## Δd - spatial separation of events
## c - speed of light

# Conditions
## - The spacetime in which the two events occur is flat (special relativity case).

# Notes
## - If `Δs**2 > 0`, the spacetime interval is said to be _timelike_. Events with a timelike separation can be causally connected,
##   i.e. one can find an inertial reference frame in which both events happen at the same place at different times.
## - If `Δs**2 = 0`, the spacetime interval is said to be _lightlike_. Events with a lightlike separation are exactly far enough
##   from each other that light could be present at both events, and they are causally connected.
## - If `Δs**2 < 0`, the spacetime interval is said to be _spacelike_. Events with a spacelike separation are causally disconnected,
##   i.e. one can find an inertial reference frame in which both events happen at the same time in different positions.

spacetime_interval = Symbol("spacetime_interval", units.length)
temporal_distance = Symbol("temporal_distance", units.time)
spatial_distance = Symbol("spatial_distance", units.length)

law = Eq(spacetime_interval**2, (speed_of_light * temporal_distance)**2 - spatial_distance**2)


@validate_input(
    temporal_distance_=temporal_distance,
    spatial_distance_=spatial_distance,
)
@validate_output(spacetime_interval)
def calculate_spacetime_interval(
    temporal_distance_: Quantity,
    spatial_distance_: Quantity,
) -> Quantity:
    expr = solve(law, spacetime_interval)[1]
    result = expr.subs({
        temporal_distance: temporal_distance_,
        spatial_distance: spatial_distance_,
    })
    return Quantity(result)
