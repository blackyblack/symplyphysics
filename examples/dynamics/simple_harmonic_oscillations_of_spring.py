#!/usr/bin/env python3

from sympy import solve, symbols, Eq
from symplyphysics import Quantity, units, print_expression, convert_to
from symplyphysics.definitions import period_from_angular_frequency as period_def
from symplyphysics.laws.dynamics import (
    period_of_spring_from_mass as spring_period_law,)
from symplyphysics.laws.kinematic import (
    displacement_in_simple_harmonic_motion as harmonic_law,
)

# Description
## A simple harmonic oscillator consists of a block of mass 2.00 kg attached to a spring of
## stiffness 100 N/m. When t = 1.00 s, the position and velocity of the block are x = 0.129 m
## and v = 3.415 m/s. (a) What is the amplitude of the oscillations? What were the (b) position
## and (c) velocity of the block at t = 0 s?

block_mass, spring_stiffness, amplitude = symbols("block_mass spring_stiffness amplitude",
    positive=True)
time, initial_time, position1, velocity1, phase_lag = symbols(
    "time initial_time position1 velocity1 phase_lag", real=True)

values = {
    block_mass: Quantity(2.00 * units.kilogram),
    spring_stiffness: Quantity(100 * units.newton / units.meter),
    initial_time: Quantity(1.00 * units.second),
    position1: Quantity(0.129 * units.meter),
    velocity1: Quantity(3.415 * units.meter / units.second),
}

# Find angular frequency of spring's oscillations

period_expr = spring_period_law.law.rhs.subs({
    spring_period_law.symbols.basic.mass: block_mass,
    spring_period_law.spring_elasticity: spring_stiffness,
})

angular_frequency_expr = solve(period_def.law,
    period_def.circular_frequency)[0].subs(period_def.period, period_expr)

# Expression of simple harmonic oscillations of spring

position_expr = harmonic_law.law.rhs.subs({
    harmonic_law.amplitude: amplitude,
    harmonic_law.angular_frequency: angular_frequency_expr,
    harmonic_law.phase_lag: phase_lag,
    harmonic_law.time: time,
})

velocity_expr = position_expr.diff(time)

# (a) Amplitude of oscillations

# Help SymPy solve the equations by using this identity
amplitude_eqn = Eq(
    position1**2 + (velocity1 / angular_frequency_expr)**2,
    position_expr**2 + (velocity_expr / angular_frequency_expr)**2,
)

amplitude_expr = solve(amplitude_eqn, amplitude)[0]
amplitude_value = convert_to(Quantity(amplitude_expr.subs(values)), units.meter).evalf(2)

print(f"Formula of amplitude of oscillations:\n{print_expression(amplitude_expr)}\n")
print(f"Amplitude of oscillations is {amplitude_value} m.\n\n")

# (b) Position and (c) velocity at zero time

# Help SymPy solve the equations by using this identity
phase_lag_eqn = Eq(
    velocity1 / position1,
    velocity_expr / position_expr,
).simplify()

phase_lag_expr = solve(phase_lag_eqn, phase_lag)[0]

print(f"Formula of phase lag:\n{print_expression(phase_lag_expr)}\n\n")

# Separate initial time and current time (time)
amplitude_expr = amplitude_expr.subs(time, initial_time)
phase_lag_expr = phase_lag_expr.subs(time, initial_time)

initial_subs = {
    amplitude: amplitude_expr,
    phase_lag: phase_lag_expr,
    time: 0,
}

initial_position_expr = position_expr.subs(initial_subs).simplify()
initial_velocity_expr = velocity_expr.subs(initial_subs).simplify()

initial_position_value = convert_to(Quantity(initial_position_expr.subs(values)),
    units.meter).evalf(2)
initial_velocity_value = convert_to(Quantity(initial_velocity_expr.subs(values)),
    units.meter / units.second).evalf(2)

print(f"Initial position of the block is {initial_position_value} m.")
print(f"Initial velocity of the block is {initial_velocity_value} m/s.")
