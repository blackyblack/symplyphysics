#!/usr/bin/env python3

from sympy import solve, symbols
from symplyphysics import (
    print_expression,
    Quantity,
    convert_to,
    units,
)
from symplyphysics.laws.kinematic import linear_velocity_from_angular_velocity_and_radius as velocity_law
from symplyphysics.definitions import period_from_angular_frequency as period_law

# Description
## The Sun is 2.3*10^4 ly (light-years) away from the center of the Milky Way galaxy. For this example
## we will approximate its movement around the galaxy's center as a circle. Its speed during its circular
## movement is 250 km/s. (a) How long does it take the Sun to make one revolution about the galactic center?
## (b) How many revolutions has the Sun completed during the last 3.5*10^8 years?

distance_to_center = symbols("distance_to_center", positive=True)
sun_speed = symbols("sun_speed", positive=True)
time_elapsed = symbols("time_elapsed", positive=True)

values = {
    distance_to_center: Quantity(2.3e4 * units.lightyear),
    sun_speed: Quantity(250.0 * units.kilometer / units.second),
    time_elapsed: Quantity(3.5e8 * units.year),
}

sun_angular_velocity = solve(
    velocity_law.law,
    velocity_law.angular_velocity,
)[0].subs({
    velocity_law.linear_velocity: sun_speed,
    velocity_law.curve_radius: distance_to_center,
})

sun_period = solve(period_law.law, period_law.period)[0].subs(
    period_law.angular_frequency,
    sun_angular_velocity,
)
sun_period_value = convert_to(
    Quantity(sun_period.subs(values)),
    units.year,
).evalf(3)

print(f"Formula for the Sun's period of revolution:\n{print_expression(sun_period)}")
print(f"The Sun makes one revolution around the Milky Way's center in {sun_period_value} years.")

count_revolutions = time_elapsed / sun_period
count_revolutions_value = Quantity(count_revolutions.subs(values)).scale_factor.evalf(3)

time_elapsed_value = convert_to(values[time_elapsed], units.year).evalf(3)

print(
    f"The Sun has made {count_revolutions_value} revolutions during the past {time_elapsed_value} years"
)
