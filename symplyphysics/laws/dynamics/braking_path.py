from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output,)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as energy_law
from symplyphysics.laws.dynamics import mechanical_work_from_force_and_move as work_law

# Description
## Let an arbitrary object move along the surface at an arbitrary speed. The friction force acts on the
## object from the surface. Then the braking path will depend on the mass of the object, its speed
## and friction force.

## Law is: S = m * v^2 / (2 * F), where
## S - braking path,
## m - mass,
## v - velocity,
## F - friction force.

braking_path = Symbol("braking_path", units.length)

mass = Symbol("mass", units.mass)
velocity = Symbol("velocity", units.velocity)
friction_force = Symbol("friction_force", units.force)

law = Eq(braking_path, mass * velocity**2 / (2 * friction_force))

# This law might be derived via "kinetic_energy_from_mass_and_velocity" law and
# "mechanical_work_from_force_and_move" law.

energy_law_applied = energy_law.law.subs({
    energy_law.body_mass: mass,
    energy_law.body_velocity: velocity
})
energy_derived = solve(energy_law_applied, energy_law.kinetic_energy_of_body, dict=True)[0][energy_law.kinetic_energy_of_body]

work_law_applied = work_law.law.subs({
    work_law.force: friction_force,
    work_law.work: energy_derived
})
distance_derived = solve(work_law_applied, work_law.distance, dict=True)[0][work_law.distance]

# Check if derived distance is same as declared.
assert expr_equals(distance_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=mass, velocity_=velocity, friction_force_=friction_force)
@validate_output(braking_path)
def calculate_braking_path(mass_: Quantity, velocity_: Quantity, friction_force_: Quantity) -> Quantity:
    result_braking_path_expr = solve(law, braking_path, dict=True)[0][braking_path]
    result_expr = result_braking_path_expr.subs({mass: mass_, velocity: velocity_, friction_force: friction_force_})
    return Quantity(result_expr)
