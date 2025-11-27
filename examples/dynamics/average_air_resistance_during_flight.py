#!/usr/bin/env python3
r"""
A body with a mass of :math:`40 \, \text{g}` was thrown vertically with an initial speed of
:math:`30 \, \frac{\text{m}}{\text{s}}`. It reached its highest point in :math:`2.5 \, \text{s}`.
Find the average force of the air's resistance exerted on the body during its flight.
"""

from sympy import solve, integrate, Idx, Eq, pi
from symplyphysics import (clone_as_function, clone_as_symbol, quantities, global_index,
    print_expression, units, convert_to_si)
from symplyphysics.definitions import (
    net_force_is_sum_of_individual_forces as superposition_law,
    acceleration_is_speed_derivative as acceleration_def,
)
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass as newtons_second_law
from symplyphysics.laws.geometry import (
    scalar_projection_is_vector_length_times_cosine_of_angle as projection_law,)
from symplyphysics.core.experimental.solvers import apply

# NOTE Since the body is moving upwards, both the force of gravity and the force of air resistance
# are directed along the negative z-axis (i.e. "down"). To simplify the solution, we project all
# vectors onto the negative z-axis.

mass = newtons_second_law.mass
time = acceleration_def.time
flight_time = clone_as_symbol(time, subscript="1")
speed = acceleration_def.speed  # projection onto -z
air_resistance = clone_as_function(newtons_second_law.force, [time])  # projection onto -z

gravity_force = solve(newtons_second_law.law, newtons_second_law.force)[0].subs({
    newtons_second_law.acceleration: quantities.acceleration_due_to_gravity,
})

net_force = superposition_law.definition.rhs.subs(global_index, Idx("i", (1, 2))).doit().subs({
    superposition_law.force[1]: gravity_force,
    superposition_law.force[2]: air_resistance(time),
})

acceleration_expr = acceleration_def.definition.rhs

newtons_second_eqn = newtons_second_law.law.subs({
    newtons_second_law.acceleration: acceleration_expr,
    newtons_second_law.force: net_force,
})

integrated_newtons_second_eqn = apply(
    newtons_second_eqn,
    lambda s: integrate(s, (time, 0, flight_time)).doit(),
)

average_air_resistance = clone_as_symbol(newtons_second_law.force, display_symbol="avg(F)")
integrated_air_resistance_expr = integrate(air_resistance(time), (time, 0, flight_time))

# TODO Should we introduce a law for averaged values?
average_air_resistance_eqn = Eq(
    average_air_resistance,
    integrated_air_resistance_expr / flight_time,
)

average_air_resistance_expr = solve(
    (integrated_newtons_second_eqn, average_air_resistance_eqn),
    (average_air_resistance, integrated_air_resistance_expr),
    dict=True,
)[0][average_air_resistance].simplify()

print(
    "Formula for the average force of air resistance:",
    print_expression(average_air_resistance_expr),
    sep="\n\n",
    end="\n\n",
)

# See NOTE above. The velocity of the body is directed along the positive z-axis.
initial_speed_qty = projection_law.law.rhs.subs({
    projection_law.vector_length: 30 * units.meter / units.second,
    projection_law.angle: pi,
})

values = {
    speed(0): initial_speed_qty,
    speed(flight_time): 0,
    mass: 40 * units.gram,
    flight_time: 2.5 * units.second,
}

average_air_resistance_qty = average_air_resistance_expr.subs(values)
average_air_resistance_si = convert_to_si(average_air_resistance_qty).evalf(2)

print(f"The average force of air resistance amounted to", average_air_resistance_si, "N.")
