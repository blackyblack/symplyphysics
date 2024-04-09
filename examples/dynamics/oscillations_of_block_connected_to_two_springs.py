#!/usr/bin/env python3

from sympy import solve, symbols
from symplyphysics import Quantity, units, convert_to, print_expression
from symplyphysics.laws.dynamics.springs import stiffness_of_two_parallel_springs as stiffness_law
from symplyphysics.laws.dynamics import period_of_spring_from_mass as spring_period_law
from symplyphysics.laws.kinematic import temporal_frequency_from_period as frequency_def

# Description
## Two identical springs of stiffness k = 7.580 kN/m are attached to a block of mass m = 0.245 kg.
## What is the frequency of oscillation on the frictionless floor?

stiffness_sym, mass_sym = symbols("stiffness, mass", positive=True)

values = {
    stiffness_sym: Quantity(7.580e3 * units.newton / units.meter),
    mass_sym: Quantity(0.245 * units.kilogram),
}

total_stiffness = solve(stiffness_law.law, stiffness_law.total_stiffness)[0].subs({
    stiffness_law.first_stiffness: stiffness_sym,
    stiffness_law.second_stiffness: stiffness_sym,
})

period = solve(spring_period_law.law, spring_period_law.oscillation_period)[0].subs({
    spring_period_law.symbols.basic.mass: mass_sym,
    spring_period_law.spring_elasticity: total_stiffness,
})

frequency = solve(frequency_def.law, frequency_def.temporal_frequency)[0].subs(
    frequency_def.period,
    period,
).simplify()

print(f"Formula of linear frequency:\n{print_expression(frequency)}\n")

frequency_quantity = Quantity(frequency.subs(values))
frequency_value = convert_to(frequency_quantity, units.hertz).evalf(3)

print(
    f"The body and two springs system oscillates with a linear frequency of {frequency_value} Hz.")
