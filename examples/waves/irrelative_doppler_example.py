#!/usr/bin/env python3

from sympy import solve
from symplyphysics import (units, convert_to, Quantity, prefixes)
from symplyphysics.laws.waves import frequency_shift_from_speed_in_collinear_motion as doppler_law

# This example show usefulness of Doppler law.
## Doppler effect is widely used in systems like radars. This effect is applicable not only to waves emitted by source, but also to waves reflected by it.
## So radar emits wave with known frequency and then receives reflected wave with frequency modified according to Doppler's effect.
## Then radar system can calculate velocity of object. Radars may use sound waves (ultrasonic), radio or infrared (electromagnetic).
## So for example ultrasonic radar emits 40 kHz wave. An object reflects this wave and at some moment radar receives this signal.
## We can calculate object's velocity.

sound_velocity = Quantity(340 * units.meter / units.second)
# Zero velocity we are going to use as a velocity of the radar related to air.
zero_velocity = Quantity(0)
emitter_frequency = Quantity(40 * prefixes.kilo * units.hertz)
# Choose any frequency and obtain the result
signal_frequency = Quantity(41.2 * prefixes.kilo * units.hertz)

solution = solve(doppler_law.law, doppler_law.source_speed,
    dict=True)[0][doppler_law.source_speed]
applied_solution = solution.subs({
    doppler_law.observed_frequency: signal_frequency,
    doppler_law.source_frequency: emitter_frequency,
    doppler_law.wave_speed: sound_velocity,
    doppler_law.observer_speed: zero_velocity
})

result_velocity = Quantity(applied_solution)
result = convert_to(result_velocity, units.kilometer / units.hour).evalf(3)

# Since object is not an emitter of signal, we cannot directly use Doppler law. Let's assume, we emit
# 40000 Hz at the approaching car. The car should observe this signal with 41200 Hz when moving at the
# resulting velocity. However, it should transmit the signal back to the radar. After applying the
# Doppler law, we should get twice the expected frequency shift - 42400 Hz.
# Therefore we should divide the resulting velocity by factor of two for the reflected signal.
result = result / 2

if result > 0:
    print(f"Object is moving away from radar with {result} km/h velocity")
elif result < 0:
    print(f"Object is moving towards radar with {-result} km/h velocity")
else:
    print("Object is not moving")
