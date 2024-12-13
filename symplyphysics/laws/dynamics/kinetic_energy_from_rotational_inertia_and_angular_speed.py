"""
Kinetic energy from rotational inertia and angular speed
========================================================

Every rotating object bears kinetic energy, which depends on the rotational inertia of the object
and its angular speed.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Kinetic_energy#Rotating_bodies>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as kinetic_energy_def
from symplyphysics.laws.kinematics import speed_via_angular_speed_and_radius as linear_velocity_law
from symplyphysics.laws.kinematics.rotational_inertia import rotational_inertia_of_particle as rotational_inertia_def

kinetic_energy = symbols.kinetic_energy
"""
The :symbols:`kinetic_energy` of the object.
"""

rotational_inertia = symbols.rotational_inertia
"""
The :symbols:`rotational_inertia` of the object.
"""

angular_speed = symbols.angular_speed
"""
The :symbols:`angular_speed` of the object.
"""

law = Eq(kinetic_energy, rotational_inertia * angular_speed**2 / 2)
"""
:laws:symbol::

:laws:latex::
"""

# Derive this law from the definition of kinetic energy and the expression for linear velocity of a rotating body
_rotation_radius = symbols.distance_to_axis

_rotational_inertia_def_subs = rotational_inertia_def.law.subs({
    rotational_inertia_def.rotational_inertia: rotational_inertia,
    rotational_inertia_def.distance_to_axis: _rotation_radius,
})
_object_mass = solve(_rotational_inertia_def_subs, rotational_inertia_def.mass)[0]

_linear_velocity_law_sub = linear_velocity_law.law.subs({
    linear_velocity_law.angular_speed: angular_speed,
    linear_velocity_law.radius_of_curvature: _rotation_radius
})
_linear_velocity = solve(_linear_velocity_law_sub, linear_velocity_law.speed)[0]

_kinetic_energy_def_sub = kinetic_energy_def.law.subs({
    kinetic_energy_def.mass: _object_mass,
    kinetic_energy_def.speed: _linear_velocity,
})
_kinetic_energy_derived = solve(_kinetic_energy_def_sub, kinetic_energy_def.kinetic_energy)[0]
_kinetic_energy_from_law = solve(law, kinetic_energy)[0]

assert expr_equals(_kinetic_energy_from_law, _kinetic_energy_derived)


@validate_input(inertia_moment_=rotational_inertia, angular_velocity_=angular_speed)
@validate_output(kinetic_energy)
def calculate_energy(inertia_moment_: Quantity, angular_velocity_: Quantity) -> Quantity:
    result_energy_expr = solve(law, kinetic_energy, dict=True)[0][kinetic_energy]
    result_expr = result_energy_expr.subs({
        rotational_inertia: inertia_moment_,
        angular_speed: angular_velocity_
    })
    return Quantity(result_expr)
