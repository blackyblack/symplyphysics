#!/usr/bin/env python3

from sympy import symbols, dsolve, solve, Eq
from symplyphysics import print_expression
from symplyphysics.laws.dynamics import (
    force_is_derivative_of_momentum as newtons_second_law,
)
from symplyphysics.laws.relativistic import (
    relativistic_momentum as relativistic_momentum_law,
)

# Description
## A particle with an initial momentum `p_0` starts moving at a moment `t = 0` due to a constant
## force `F`. Find the time dependence of the particle's velocity and the distance covered.

# Momentum is expressible as a function of time

force, initial_momentum = symbols("force, initial_momentum", positive=True)

time = newtons_second_law.time
momentum_function = newtons_second_law.momentum

momentum_via_time = dsolve(
    newtons_second_law.law.subs(newtons_second_law.force(time), force),
    momentum_function(time),
    ics={momentum_function(0): initial_momentum},
).rhs

# We can find velocity as a function of momentum

speed, rest_mass, momentum = symbols("speed, rest_mass, momentum", positive=True)

relativistic_momentum_eqn = relativistic_momentum_law.law.subs({
    relativistic_momentum_law.momentum: momentum,
    relativistic_momentum_law.velocity: speed,
    relativistic_momentum_law.rest_mass: rest_mass,
})

speed_via_momentum = solve(relativistic_momentum_eqn, speed)[0]

# Now we can express speed as a function of time

speed_via_time = speed_via_momentum.subs(momentum, momentum_via_time)

print("Formula of particle speed as a function of time:\n")
print(print_expression(Eq(speed, speed_via_time)))

# Integrate the velocity formula over time to obtain the distance formula.
# Note that the distance covered is zero at `t = 0`.

distance, length_term = symbols("distance, length_term", positive=True)

distance_via_time = speed_via_time.integrate(time)

length_term_expr = distance_via_time.subs(time, 0)

distance_via_time = solve(
    (
        Eq(distance, distance_via_time - length_term_expr),
        Eq(length_term, length_term_expr),
    ),
    (distance, rest_mass),
    dict=True,
)[0][distance]

print("\n\nFormula of particle distance as a function of time:")
print(print_expression(Eq(distance, distance_via_time)))

print("\n\nFormula of `length_term`:")
print(print_expression(length_term_expr))

# Now let us look at the edge case `F = 0`. That is, the particle continues its movement as if no force has ever
# acted upon it.

zero_force_momentum = momentum_via_time.subs(force, 0)

zero_force_speed = speed_via_momentum.subs(momentum, zero_force_momentum)

print("\n\nFormula of particle speed when force is zero:\n")
print(print_expression(zero_force_speed))

zero_force_distance = zero_force_speed.integrate(time)

print("\n\nFormula of particle distance when force is zero:\n")
print(print_expression(zero_force_distance))
