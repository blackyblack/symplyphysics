from sympy import Eq, sin, pi
from symplyphysics import (
    units,
    angle_type,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)
from symplyphysics.core.symbols.quantities import scale_factor

# Description
## Suppose a body is falling freely in Earth's gravity field with its initial velocity being zero.
## Then the effect of the Coriolis force on the falling body can be found in the fact that it deflects
## from plumbline in the easterly and southerly (equatorial) directions.

# Conditions
## - The vector of free fall acceleration `g` is considered constant.

# Notes
## - The southerly deviation is extremely small and almost unobservable due to the `t / T` factor.

# Law: s_south = pi * (t / T) * s_east * sin(theta)
## s_sourth - southerly (equatorial) deviation of falling body from plumbline due to Earth's rotation
## t - time of the body's fall
## T - period of Earth's rotation
## s_east - easterly deviation of falling body from plumbline
## theta - latitude of the place the body is located in

southerly_deviation_from_plumbline = Symbol("southerly_deviation_from_plumbline", units.length)
fall_time = Symbol("fall_time", units.time)
rotation_period = Symbol("rotation_period", units.time)
easterly_deviation_from_plumbline = Symbol("easterly_deviation_from_plumbline", units.length)
latitude = Symbol("latitude", angle_type)

law = Eq(
    easterly_deviation_from_plumbline,
    pi * (fall_time / rotation_period) * easterly_deviation_from_plumbline * sin(latitude),
)

# TODO: derive from the solution of the relative motion equation in Earth's gravitational field.


@validate_input(
    fall_time_=fall_time,
    rotation_period_=rotation_period,
    easterly_deviation_from_plumbline_=easterly_deviation_from_plumbline,
    latitude_=latitude,
)
@validate_output(easterly_deviation_from_plumbline)
def calculate_southerly_deviation_from_plumbline(
    fall_time_: Quantity,
    rotation_period_: Quantity,
    easterly_deviation_from_plumbline_: Quantity,
    latitude_: Quantity | float,
) -> Quantity:
    result = law.rhs.subs({
        fall_time: fall_time_,
        rotation_period: rotation_period_,
        easterly_deviation_from_plumbline: easterly_deviation_from_plumbline_,
        latitude: scale_factor(latitude_),
    })
    return Quantity(result)
