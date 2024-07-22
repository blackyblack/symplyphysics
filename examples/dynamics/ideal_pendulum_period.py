#!/usr/bin/env python3

from symplyphysics import Quantity, convert_to, units, print_expression
from symplyphysics.laws.dynamics import period_of_ideal_pendulum_from_length as pendulum_period
from symplyphysics.definitions import temporal_frequency_from_period as frequency_def

# This example calculates ideal pendulum period from its length

print(f"Pendulum oscillation period formula is:\n{print_expression(pendulum_period.law)}")

pendulum_length = Quantity(24 * units.inch)

oscillating_period_expr = pendulum_period.law.subs(pendulum_period.pendulum_length,
    pendulum_length).rhs
oscillating_period = Quantity(oscillating_period_expr)
oscillating_period_seconds = convert_to(oscillating_period, units.second).evalf(5)

pendulum_length_inches = int(convert_to(pendulum_length, units.inch).evalf(3))
print(
    f"Pendulum with length = {pendulum_length_inches} {units.inch} has oscillating period {oscillating_period_seconds} {units.second}\n"
)

oscillating_frequency_hertz = frequency_def.law.rhs.subs(frequency_def.period,
    oscillating_period_seconds)
print(f"Oscillating frequency is {oscillating_frequency_hertz} {units.hertz}\n")
