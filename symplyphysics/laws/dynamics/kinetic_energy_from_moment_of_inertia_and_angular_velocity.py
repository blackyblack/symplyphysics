from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity, SI
)
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

# Description
## If an object has a inertia moment and spins with some angular velocity, it bears kinetic energy.
## Law: E = I * w**2 / 2, where
## E is kinetic energy of spinning object
## I is inertia moment of this object
## w is angular velocity

kinetic_energy, object_inertia_moment, angular_velocity = symbols('kinetic_energy object_inertia_moment angular_velocity')
law = Eq(kinetic_energy, object_inertia_moment * angular_velocity**2 / 2)


def print():
    return pretty(law, use_unicode=False)

@validate_input(inertia_moment_=units.mass * units.length**2, angular_velocity_=angle_type / units.time)
@validate_output(units.energy)
def calculate_energy(inertia_moment_: Quantity, angular_velocity_: Quantity) -> Quantity:
    result_energy_expr = solve(law, kinetic_energy, dict=True)[0][kinetic_energy]
    #HACK: sympy angles are always in radians and angular velocity cannot be properly converted to velocity
    SI.set_quantity_dimension(angular_velocity_, 1 / units.time)
    result_expr = result_energy_expr.subs({object_inertia_moment: inertia_moment_, angular_velocity: angular_velocity_})
    return expr_to_quantity(result_expr, 'kinetic_energy')

