"""
Thin rod about axis through center perpendicular to length
==========================================================

A thin rod is rotating about an axis that passes through its center perpendicular to its length.

**Conditions:**

#. The rod is uniform.
#. The rod is thin, i.e. its diameter is much less than its length.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics.rotational_inertia.geometries import (
    slab_about_perpendicular_axis_through_center as slab_formula)

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
"""
Rotational inertia of the rod.

Symbol:
    :code:`I`
"""

length = Symbol("length", units.length)
"""
Length of the rod.

Symbol:
    :code:`l`
"""

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of the rod.

Symbol:
    :code:`m`
"""

law = Eq(rotational_inertia, mass * length**2 / 12)
r"""
:code:`I = 1/12 * m * l^2`

Latex:
    .. math::
        I = \frac{1}{12} m l^2
"""

# Derive law from formula for a slab rotating about the axis perpendicular to its length and width
# passing through its center. The thin rod is a particular case of it, when the width of the slab
# approaches zero.

_rotational_inertia_derived = slab_formula.law.rhs.subs({
    slab_formula.mass: mass,
    slab_formula.length: length,
    slab_formula.width: 0,
})

assert expr_equals(law.rhs, _rotational_inertia_derived)


@validate_input(mass_=mass, length_=length)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(mass_: Quantity, length_: Quantity) -> Quantity:
    result = law.rhs.subs({
        mass: mass_,
        length: length_,
    })
    return Quantity(result)
