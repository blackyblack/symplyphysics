from sympy import Eq, cos, pi
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

# Law: s_east = (4 * pi / 3) * (t / T) * h * cos(theta)
## s_east - easterly deviation of falling body from plumbline due to Earth's rotation
## t - time of the body's fall
## T - period of Earth's rotation
## h - initial elevation of the body from Earth's surface
## theta - latitude of the place the body is located in

easterly_deviation_from_plumbline = Symbol("easterly_deviation_from_plumbline", units.length)
fall_time = Symbol("fall_time", units.time)
rotation_period = Symbol("rotation_period", units.time)
initial_elevation = Symbol("initial_elevation", units.length)
latitude = Symbol("latitude", angle_type)

law = Eq(
    easterly_deviation_from_plumbline,
    (4 * pi / 3) * (fall_time / rotation_period) * initial_elevation * cos(latitude),
)

# TODO: derive from the solution of the relative motion equation in Earth's gravitational field.


@validate_input(
    fall_time_=fall_time,
    rotation_period_=rotation_period,
    initial_elevation_=initial_elevation,
    latitude_=latitude,
)
@validate_output(easterly_deviation_from_plumbline)
def calculate_easterly_deviation_from_plumbline(
    fall_time_: Quantity,
    rotation_period_: Quantity,
    initial_elevation_: Quantity,
    latitude_: Quantity | float,
) -> Quantity:
    result = law.rhs.subs({
        fall_time: fall_time_,
        rotation_period: rotation_period_,
        initial_elevation: initial_elevation_,
        latitude: scale_factor(latitude_),
    })
    return Quantity(result)
