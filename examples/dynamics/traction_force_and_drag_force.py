from sympy import Idx, solve, Symbol, Eq
from symplyphysics import print_expression, global_index
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass as second_newton_law
from symplyphysics.laws.kinematics import position_via_constant_acceleration_and_time as distance_law
from symplyphysics.definitions import net_force_is_sum_of_individual_forces as superposition_law

# A trolleybus with a mass of 12 tons, starting from a place,
# passes a distance of 10 m along a horizontal path in 5 seconds.
# Determine the thrust force developed by the engine if the drag force is 2.4 kN.

time_of_motion = Symbol("time_of_motion")
mass_of_trolleybus = Symbol("mass_of_trolleybus")
drag_force = Symbol("drag_force")
distance = Symbol("distance")

traction_force = Symbol("traction_force")

# The traction force is directed horizontally in the direction of movement
# of the trolleybus, and the drag force is directed horizontally against
# the movement.
index_local = Idx("index_local", (1, 2))
superposition_of_two_forces = superposition_law.definition.subs(global_index, index_local).doit()

acceleration_force = superposition_of_two_forces.subs({
    superposition_law.force[1]: traction_force,
    superposition_law.force[2]: -drag_force,
}).rhs

acceleration_value = second_newton_law.law.subs({
    second_newton_law.mass: mass_of_trolleybus,
    second_newton_law.force: acceleration_force
}).rhs

distance_value = distance_law.law.subs({
    distance_law.time: time_of_motion,
    distance_law.acceleration: acceleration_value,
    distance_law.initial_speed: 0,
    distance_law.initial_position: 0,
    distance_law.final_position: distance
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
