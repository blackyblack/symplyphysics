"""
Slab about perpendicular axis through center
============================================

A solid slab is rotating about an axis that is perpendicular to the slab and goes through its center.

**Conditions:**

#. The slab is uniform.
"""

from sympy import Eq, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics.rotational_inertia import (
    rotational_inertia_cartesian_integral as integral_law,)
from symplyphysics.definitions import density_from_mass_volume as density_def

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
"""
Rotational inertia of the slab.

Symbol:
    :code:`I`
"""

length = Symbol("length", units.length)
"""
Length, or first slab dimension perpendicular to the axis.

Symbol:
    :code:`a`
"""

width = Symbol("width", units.length)
"""
Width, or second slab dimension perpendicular to the axis.

Symbol:
    :code:`b`
"""

mass = symbols.mass
"""
The :attr:`~symplyphysics.symbols.mass` of the slab.
"""

law = Eq(rotational_inertia, mass * (length**2 + width**2) / 12)
r"""
:code:`I = 1/12 * m * (a^2 + b^2)`

Latex:
    .. math::
        I = \frac{1}{12} m (a^2 + b^2)
"""

# Derive this law from the integral definition of rotational inertia in cartesian coordinates.
# Condition: _density of slab is constant.

# Reference frame:
## z-axis is parallel to the rotational axis in question (_height of slab)
## x-axis and y-axis are perpendicular to the rotational axis (length and width of slab)

_height = Symbol("_height", units.length)
_volume = length * width * _height

_density = density_def.definition.rhs.subs({
    density_def.mass: mass,
    density_def.volume: _volume,
})

_distance_to_axis = sqrt(integral_law.x**2 + integral_law.y**2)

_rotational_inertia_derived = integral_law.law.rhs.subs({
    integral_law.density(integral_law.x, integral_law.y, integral_law.z):
        _density,
    integral_law.distance_to_axis(integral_law.x, integral_law.y, integral_law.z):
        _distance_to_axis,
    integral_law.x_start:
    -1 * length / 2,
    integral_law.y_start:
    -1 * width / 2,
    integral_law.z_start:
        0,
    integral_law.x_end:
    length / 2,
    integral_law.y_end:
    width / 2,
    integral_law.z_end:
        _height,
}).doit().simplify()

assert expr_equals(law.rhs, _rotational_inertia_derived)


@validate_input(mass_=mass, length_=length, width_=width)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, length_: Quantity, width_: Quantity) -> Quantity:
    result = law.rhs.subs({
        mass: mass_,
        length: length_,
        width: width_,
    })
    return Quantity(result)
