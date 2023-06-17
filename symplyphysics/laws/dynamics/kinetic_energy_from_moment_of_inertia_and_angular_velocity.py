from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression, angle_type,
    validate_input, validate_output)

# Description
## If an object has a inertia moment and spins with some angular velocity, it bears kinetic energy.
## Law: E = I * w**2 / 2, where
## E is kinetic energy of spinning object
## I is inertia moment of this object
## w is angular velocity

kinetic_energy = Symbol("kinetic_energy", units.energy)
object_inertia_moment = Symbol("object_inertia_moment", units.mass * units.length**2)
angular_velocity = Symbol("angular_velocity", angle_type / units.time)

law = Eq(kinetic_energy, object_inertia_moment * angular_velocity**2 / 2)


def print() -> str:
    return print_expression(law)


@validate_input(inertia_moment_=object_inertia_moment, angular_velocity_=angular_velocity)
@validate_output(kinetic_energy)
def calculate_energy(inertia_moment_: Quantity, angular_velocity_: Quantity) -> Quantity:
    result_energy_expr = solve(law, kinetic_energy, dict=True)[0][kinetic_energy]
    result_expr = result_energy_expr.subs({
        object_inertia_moment: inertia_moment_,
        angular_velocity: angular_velocity_
    })
    return expr_to_quantity(result_expr)
