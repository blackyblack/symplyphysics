"""
Solid disk about central axis
=============================

A solid disk (cylinder) rotates about its central axis (axis of cylindrical symmetry).

**Conditions:**

#. The disk is uniform.
"""

from sympy import Eq, pi
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics.rotational_inertia import rotational_inertia_cylindrical_integral as integral_law
from symplyphysics.definitions import density_from_mass_volume as density_def

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
"""
Rotational inertia of the disk.

Symbol:
    :code:`I`
"""

radius = Symbol("radius", units.length)
"""
Radius of the disk.

Symbol:
    :code:`r`
"""

mass = symbols.mass
"""
The :symbols:`mass` of the disk.
"""

law = Eq(rotational_inertia, mass * radius**2 / 2)
r"""
:code:`I = 1/2 * m * r^2`

Latex:
    .. math::
        I = \frac{1}{2} m r^2
"""

# Derive law from general integral in cylindrical coordinates

_length = Symbol("_length", units.length)
_volume = pi * radius**2 * _length

_density = density_def.definition.rhs.subs({
    density_def.mass: mass,
    density_def.volume: _volume,
})

_density_applied_sym = integral_law.density(integral_law.radius, integral_law.polar_angle,
    integral_law.height)

_rotational_inertia_derived = integral_law.law.rhs.subs({
    _density_applied_sym: _density,
    integral_law.radius_start: 0,
    integral_law.radius_end: radius,
    integral_law.polar_angle_start: 0,
    integral_law.polar_angle_end: 2 * pi,
    integral_law.height_start: 0,
    integral_law.height_end: _length,
}).doit()

assert expr_equals(_rotational_inertia_derived, law.rhs)


@validate_input(mass_=mass, radius_=radius)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, radius_: Quantity) -> Quantity:
    result = law.rhs.subs({
        mass: mass_,
        radius: radius_,
    })
    return Quantity(result)
