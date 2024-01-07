from sympy import solve, Symbol, Eq
from symplyphysics import print_expression
from symplyphysics.laws.dynamics import acceleration_from_force as second_newton_law
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as distance_law

# A trolleybus with a mass of 12 tons, starting from a place,
# passes a distance of 10 m along a horizontal path in 5 seconds.
# Determine the thrust force developed by the engine if the drag force is 2.4 kN.

time_of_motion = Symbol("time_of_motion")
mass_of_trolleybus = Symbol("mass_of_trolleybus")
drag_force = Symbol("drag_force")
distance = Symbol("distance")

traction_force = Symbol("traction_force")

# According to Newton's second law,
# the thrust force is directed horizontally in the direction of movement
# of the trolleybus, and the drag force is directed horizontally against
# the movement. In the projection on the horizontal axis,
# the sum of the forces acting on the trolleybus
# F = F_trac - F_drag
acceleration_value = second_newton_law.law.subs({
    second_newton_law.mass: mass_of_trolleybus,
    second_newton_law.force: traction_force - drag_force
}).rhs

distance_value = distance_law.law.subs({
    distance_law.movement_time: time_of_motion,
    distance_law.constant_acceleration: acceleration_value,
    distance_law.initial_velocity: 0,
    distance_law.distance(distance_law.movement_time): distance
})
print(f"Final equation:\n{print_expression(distance_value)}")
traction_force_equation = solve(distance_value, traction_force, dict=True)[0][traction_force]
answer = Eq(traction_force, traction_force_equation)
print(f"Total traction force equation:\n{print_expression(answer)}")
traction_force_N = traction_force_equation.subs({
    time_of_motion: 5,
    distance: 10,
    mass_of_trolleybus: 12_000,
    drag_force: 2_400
})
print(f"Traction force is: {traction_force_N} N")
