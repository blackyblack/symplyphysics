from collections import namedtuple

from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity, SI
)
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

# Description
## If an object has a inertia moment and spins with some angular velocity, it bears kinetic energy.
## Law: E = I * w**2 / 2, where
## E is kinetic energy of moving object
## I is inertia moment of this object
## V is angular velocity

kinetic_energy, object_inertia_moment, angular_velocity = symbols('kinetic_energy object_inertia_moment angular_velocity')
law = Eq(kinetic_energy, object_inertia_moment * angular_velocity**2 / 2)

I = units.Quantity('I')
SI.set_quantity_dimension(I, units.mass * units.length**2)
SI.set_quantity_scale_factor(I, 2 * units.kilogram * units.meter**2)

w = units.Quantity('w')
SI.set_quantity_dimension(w, angle_type / units.time)
SI.set_quantity_scale_factor(w, 3 * units.radians / units.second)

ang_velocity_radians = w.scale_factor
result_energy_expr = solve(law, kinetic_energy, dict=True)[0][kinetic_energy]
print(f"{result_energy_expr}")
result_expr = result_energy_expr.subs({object_inertia_moment: I, angular_velocity: ang_velocity_radians})
print(f"{result_expr}")

def print():
    return pretty(law, use_unicode=False)

@validate_input(moment_ = units.mass * units.length**2, ang_velocity_= angle_type / units.time)
@validate_output(units.energy)
def calculate_energy(moment_: Quantity, ang_velocity_: Quantity) -> Quantity:
    ang_velocity_radians = ang_velocity_.scale_factor    
    result_energy_expr = solve(law, kinetic_energy, dict=True)[0][kinetic_energy]
    result_expr = result_energy_expr.subs({object_inertia_moment: moment_, angular_velocity: ang_velocity_radians})
    return expr_to_quantity(result_expr, 'energy')

