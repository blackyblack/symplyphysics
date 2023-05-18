#!/usr/bin/env python3

from symplyphysics import print_expression
from sympy import (Derivative, Function, cos, diff, sin, sqrt, symbols, Eq, simplify)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector
from symplyphysics.laws.dynamics import potential_energy_from_mass_and_height as potential_energy
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy
from symplyphysics.laws.kinematic import linear_velocity_from_angular_velocity_and_radius as angular_velocity_law
from symplyphysics.laws.kinematic import period_from_circular_frequency as angular_frequency

#TODO: move it to laws.
#TODO: add some real numbers to make it a proper example

# This example calculates ideal pendulum period from its length, mass and free fall acceleration.
## Ideal pendulum is an object hanging on a thread. In a field of gravitation it starts oscillating after been pushed out of balance.

# Conditions:
## 1. Ideal pendulum doesn't accept or loose any energy. No any friction.
## 2. Object is small.
## 3. Another end of a thread is not moving in current coordinate system.
## 4. Thread is weightless and doesn't change its length.

pendulum_length = symbols("pendulum_length")
pendulum_mass = symbols("pendulum_mass")
time = symbols("time")

## 2-dimension system is selected for this task. Y-axis is along gravity vector. X-axis directed right. Zero is in balanced position.

## Pendulum angle is angle between thread and gravity vector. In balanced position it is 0.
pendulum_angle = symbols("pendulum_angle", cls=Function)

## There are two forces in this task: thread reaction force as centripetal force and gravity force. Gravity force causes free fall acceleration - g independently from mass.
## Projection of gravity force to a tangental velocity causes tangential acceleration.
gravity_acceleration = symbols("gravity_acceleration")

## Pendulum oscillation is cyclic transfer of energy from kinetic to potential. To set oscillation up we have to input some energy. Usually it is done by biasing the pendulum to some angle and letting it go.
## Biasing the pendulum is giving to it some amount of potential energy.

pendulum_height_before = pendulum_length
pendulum_height_after = projector.law.subs({
    projector.vector_length: pendulum_length,
    projector.vector_angle: pendulum_angle(time)
}).rhs
amount_of_potential_energy = potential_energy.law.subs({
    potential_energy.body_mass: pendulum_mass,
    potential_energy.free_fall_acceleration: gravity_acceleration,
    potential_energy.height: (pendulum_height_before - pendulum_height_after)
}).rhs

## Kinetic energy of the pendulum is:
## pendulum_mass * (pendulum_length * angular_velocity)**2 / 2

linear_velocity = angular_velocity_law.law.subs({
    angular_velocity_law.curve_radius: pendulum_length,
    angular_velocity_law.angular_velocity: Derivative(pendulum_angle(time), time),
}).rhs
amount_of_kinetic_energy = kinetic_energy.law.subs({
    kinetic_energy.body_mass: pendulum_mass,
    kinetic_energy.body_velocity: linear_velocity
}).rhs

# Total energy is constant
total_energy = symbols("total_energy", constant=True)
time = symbols("time")

total_energy_eq = Eq(total_energy, amount_of_kinetic_energy + amount_of_potential_energy)

# Differentiate both sides of equation
total_energy_diff_eq = Eq(diff(total_energy_eq.lhs, time), diff(total_energy_eq.rhs, time))
denominator = pendulum_mass * pendulum_length * Derivative(pendulum_angle(time), time)

# Will result in harmonic oscillator equation:
## Derivative(pendulum_angle(time), (time, 2)) + gravity_acceleration / pendulum_length * sin(pendulum_angle(time)) = 0
harmonic_oscillation_eq = Eq(simplify(total_energy_diff_eq.lhs / denominator), simplify(total_energy_diff_eq.rhs / denominator))

#NOTE: large displacement angle (over 15 degrees) gives quite a complex solution for the differential equation.

# For small angles, sin(pendulum_angle) can be reduced to pendulum_angle
small_angle_harmonic_oscillation = harmonic_oscillation_eq.subs(sin(pendulum_angle(time)), pendulum_angle(time))

# dsolve() gives us solution in exponential form - we are looking for the solution as trigonometric function
maximum_angle = symbols("maximum_angle", constant=True)
initial_angle = symbols("initial_angle", constant=True)
oscillation_angular_frequency = sqrt(gravity_acceleration / pendulum_length)
angle_function_solution = maximum_angle * cos(oscillation_angular_frequency * time + initial_angle)
dsolved = small_angle_harmonic_oscillation.subs(pendulum_angle(time), angle_function_solution)
assert expr_equals(dsolved.lhs, dsolved.rhs)

# angle_function_solution is a periodic function with angular frequency equal to oscillation_angular_frequency

print("Pendulum oscillation angular frequency formula is:\n{}".format(print_expression(oscillation_angular_frequency)))

oscillation_period = angular_frequency.law.subs(angular_frequency.circular_frequency, oscillation_angular_frequency)

print("Pendulum oscillation period formula is:\n{}".format(print_expression(oscillation_period)))
