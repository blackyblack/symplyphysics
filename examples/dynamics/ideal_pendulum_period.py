#!/usr/bin/env python3

from sympy import solve, symbols, sin, cos, pi, Eq
from symplyphysics import (print_expression, units, convert_to, Quantity)
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector

# This example calculates ideal pendulum period from it's length, mass and free fall acceleration.
## Ideal pendulum is an object hanging on a thread. In a field of gravitation it starts oscillating after been pushed out of balance.

# Conditions:
## 1. Ideal pendulum doesnt accept or loose any energy. No any friction.
## 2. Object is small.
## 3. Another end of a thread is not moving in current axis.
## 4. Thread is weightless and doesnt change it's length.

pendulum_length = symbols("pendulum_length")
pendulum_mass = symbols("pendulum_mass")

## 2-dimension system is selected for this task. Y-axis is along gravity vector. X-axis directed right. Zero is in balanced position.

## Pendulum angle is angle between thread and gravity vector. In balanced position it is 0.
pendulum_angle = symbols("pendulum_angle")

## There are two forces in this task: thread reaction force as centripetal force and gravity force. Gravity force causes free fall acceleration - g independently from mass.
## Projection of gravity force to a tangental velocity causes tangential acceleration.
gravity_acceleration = symbols("gravity_acceleration")
tangential_acceleration = symbols("tangential_acceleration")

tangential_acceleration = gravity_acceleration * sin(pi/2 - pendulum_angle )

## Pendulum oscillation is cyclic transfer of energy from kinetic to potential. To set oscillation up we have to input some energy. Usually it is done by biasing the pendulum to some angle and letting it go.
## Biasing the pendulum is giving to it some amount of potential energy.

magnitude_angle = symbols("magnitude_angle")

amount_of_potential_energy = symbols("amount_of_potential_energy")

amount_of_potential_energy = pendulum_mass * gravity_acceleration * (pendulum_length - pendulum_length * cos(magnitude_angle))

## After 1/4 of oscillating period pendulum reaches it's lowest position with angle = 0 and maximum of velocity (all potential energy has been transformed to kinetic)

magnitude_velocity = symbols("magnitude_velocity")
amount_of_kinetic_energy = pendulum_mass * magnitude_velocity * magnitude_velocity / 2

Equation_1 = Eq(amount_of_kinetic_energy, amount_of_potential_energy)

## The pendulum's movement in 1/4 period is arc as a proportional part of 2*pi arc

path = 2 * pi * pendulum_length * pendulum_angle / (2 * pi)

## Also path is double-integrated acceleration with zero initial velocity.

Equation_2 = Eq(path, integral(integral(tangential_acceleration)))

## So that's the question. What is integration variable - time or angle? And what next?