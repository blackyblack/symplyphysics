#!/usr/bin/env python3

from symplyphysics import Quantity, convert_to, expr_to_quantity, print_expression, units
from symplyphysics.laws.dynamics import period_of_ideal_pendulum_from_length as pendulum_period

# This example calculates ideal pendulum period from its length

print("Pendulum oscillation period formula is:\n{}".format(print_expression(pendulum_period.law)))

pendulum_length = Quantity(24 * units.inch)
pendulum_length_inches = int(convert_to(pendulum_length, units.inch).subs(units.inch, 1).evalf(3))

oscillating_period_expr = pendulum_period.law.subs(pendulum_period.pendulum_length,
    pendulum_length).rhs
oscillating_period = expr_to_quantity(oscillating_period_expr)
oscillating_period_seconds = convert_to(oscillating_period, units.second).subs(units.second,
    1).evalf(5)

print(
    f"Pendulum with length = {pendulum_length_inches} {units.inch} has oscillating period {oscillating_period_seconds} {units.second}\n"
)

#TODO: add frequency definition

oscillating_frequency_hertz = 1 / oscillating_period_seconds

print(f"Oscillating frequency is {oscillating_frequency_hertz} {units.hertz}\n")
