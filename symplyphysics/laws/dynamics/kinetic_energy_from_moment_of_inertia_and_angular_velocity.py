from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    angle_type,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy_def
from symplyphysics.laws.kinematic import linear_velocity_from_angular_velocity_and_radius as linear_velocity_law
from symplyphysics.laws.kinematic.rotational_inertia import rotational_inertia_of_particle as rotational_inertia_def

# Description
## If an object has a inertia moment and spins with some angular velocity, it bears kinetic energy.
## Law: E = I * w**2 / 2, where
## E is kinetic energy of spinning object
## I is inertia moment of this object
## w is angular velocity

kinetic_energy = Symbol("kinetic_energy", units.energy)
object_inertia_moment = Symbol("object_inertia_moment", units.mass * units.area)
angular_velocity = Symbol("angular_velocity", angle_type / units.time)

law = Eq(kinetic_energy, object_inertia_moment * angular_velocity**2 / 2)

# Derive this law from the definition of kinetic energy and the expression for linear velocity of a rotating body
rotation_radius = Symbol("rotation_radius", units.length)

rotational_inertia_def_subs = rotational_inertia_def.law.subs({
    rotational_inertia_def.rotational_inertia: object_inertia_moment,
    rotational_inertia_def.radius: rotation_radius,
})
object_mass = solve(rotational_inertia_def_subs, rotational_inertia_def.particle_mass)[0]

linear_velocity_law_sub = linear_velocity_law.law.subs({
    linear_velocity_law.angular_velocity: angular_velocity,
    linear_velocity_law.curve_radius: rotation_radius
})
linear_velocity = solve(linear_velocity_law_sub, linear_velocity_law.linear_velocity)[0]

kinetic_energy_def_sub = kinetic_energy_def.law.subs({
    kinetic_energy_def.symbols.basic.mass: object_mass,
    kinetic_energy_def.body_velocity: linear_velocity,
})
kinetic_energy_derived = solve(kinetic_energy_def_sub, kinetic_energy_def.kinetic_energy_of_body)[0]
kinetic_energy_from_law = solve(law, kinetic_energy)[0]

assert expr_equals(kinetic_energy_from_law, kinetic_energy_derived)


def print_law() -> str:
    return print_expression(law)


@validate_input(inertia_moment_=object_inertia_moment, angular_velocity_=angular_velocity)
@validate_output(kinetic_energy)
def calculate_energy(inertia_moment_: Quantity, angular_velocity_: Quantity) -> Quantity:
    result_energy_expr = solve(law, kinetic_energy, dict=True)[0][kinetic_energy]
    result_expr = result_energy_expr.subs({
        object_inertia_moment: inertia_moment_,
        angular_velocity: angular_velocity_
    })
    return Quantity(result_expr)
