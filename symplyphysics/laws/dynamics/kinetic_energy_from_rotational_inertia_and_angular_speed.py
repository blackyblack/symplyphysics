"""
Kinetic energy from rotational inertia and angular speed
========================================================

Every rotating object bears kinetic energy, which depends on the rotational inertia of the object
and its angular speed.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    angle_type,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as kinetic_energy_def
from symplyphysics.laws.kinematic import linear_velocity_from_angular_velocity_and_radius as linear_velocity_law
from symplyphysics.laws.kinematic.rotational_inertia import rotational_inertia_of_particle as rotational_inertia_def

kinetic_energy = Symbol("kinetic_energy", units.energy)
r"""
The kinetic energy of the object.

Symbol:
    :code:`K`
"""

rotational_inertia = Symbol("rotational_inertia", units.mass * units.area)
"""
The rotational inertia of the object.

Symbol:
    :code:`I`
"""

angular_speed = Symbol("angular_speed", angle_type / units.time)
r"""
The angular speed of the object.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

law = Eq(kinetic_energy, rotational_inertia * angular_speed**2 / 2)
r"""
:code:`K = I * w^2 / 2`

Latex:
    .. math::
        K = \frac{1}{2} I \omega^2
"""

# Derive this law from the definition of kinetic energy and the expression for linear velocity of a rotating body
_rotation_radius = Symbol("_rotation_radius", units.length)

_rotational_inertia_def_subs = rotational_inertia_def.law.subs({
    rotational_inertia_def.rotational_inertia: rotational_inertia,
    rotational_inertia_def.radius: _rotation_radius,
})
_object_mass = solve(_rotational_inertia_def_subs, rotational_inertia_def.particle_mass)[0]

_linear_velocity_law_sub = linear_velocity_law.law.subs({
    linear_velocity_law.angular_velocity: angular_speed,
    linear_velocity_law.curve_radius: _rotation_radius
})
_linear_velocity = solve(_linear_velocity_law_sub, linear_velocity_law.linear_velocity)[0]

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
