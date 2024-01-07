from sympy import solve, Symbol, Eq
from symplyphysics import print_expression
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projection_velocity
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy

# The initial velocity of the bullet is 600 m/s, its mass is 10 g.
# At what angle to the horizon did it fly out of the muzzle of the gun,
# if its kinetic energy at the highest point of the trajectory is 450 J?

start_velocity = Symbol("start_velocity")
kinetic_energy_in_pick = Symbol("kinetic_energy_in_pick")
mass_of_bullet = Symbol("mass_of_bullet")

angle_of_shot = Symbol("angle_of_shot")

velocity_projection_equation = projection_velocity.law.subs(({
    projection_velocity.vector_length: start_velocity,
    projection_velocity.vector_angle: angle_of_shot
})).rhs

kinetic_energy_equation = kinetic_energy.law.subs({
    kinetic_energy.kinetic_energy_of_body: kinetic_energy_in_pick,
    kinetic_energy.body_velocity: velocity_projection_equation,
    kinetic_energy.body_mass: mass_of_bullet
})
print(f"Final equation: {print_expression(kinetic_energy_equation)}")

angle_of_shot_equation = solve(kinetic_energy_equation, angle_of_shot, dict=True)[0][angle_of_shot]
answer = Eq(angle_of_shot, angle_of_shot_equation)
print(f"Total angle of shot equation:\n{print_expression(answer)}")

angle_of_shot_rad = angle_of_shot_equation.subs({
    start_velocity: 600,
    mass_of_bullet: 0.01,
    kinetic_energy_in_pick: 450
})
print(f"Angle of shot is: {angle_of_shot_rad} rad")
