"""
Magnetic field due to infinite wire
===================================

The magnitude of the magnetic flux density due to a thin, straight, infinite wire depends on the
current through it and the radial distance to the wire.

**Notation:**

#. :quantity_notation:`vacuum_permeability`.

**Conditions:**

#. The wire is uniform, straight, and thin.

#. The vector of the magnetic flux density is oriented in space according to the right-hand rule.

**Links:**

#. `Physics LibreTexts, formula in the box <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/12%3A_Sources_of_Magnetic_Fields/12.03%3A_Magnetic_Field_due_to_a_Thin_Straight_Wire>`__.
"""

from sympy import (Eq, pi, atan2, S, limit)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)
from symplyphysics.laws.electricity import magnetic_induction_of_linear_conductor_of_finite_length as _finite_wire_law

from symplyphysics.core.expr_comparisons import expr_equals

magnetic_flux_density = symbols.magnetic_flux_density
"""
Magnitude of :symbols:`magnetic_flux_density`.
"""

absolute_permeability = symbols.absolute_permeability
"""
:symbols:`absolute_permeability` of the medium around the wire.
"""

current = symbols.current
"""
:symbols:`current` flowing through the wire.
"""

radial_distance = symbols.distance_to_axis
"""
Radial distance to wire. See :symbols:`distance_to_axis`.
"""

law = Eq(
    magnetic_flux_density,
    absolute_permeability * current / (2 * pi * radial_distance),
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive this law from the case of a finite conductor.

# When the length of a finite conductor approaches infinity, the angle made between the wire and
# the vector directed from the point of interest to a wire's end approaches zero from both sides
# of the wire. Let `R` is the radial distance to the wire and `x` the distance to the wire's end
# measured along the wire, then the angle `theta` is given by:

_r = symbols.distance_to_axis
_x = symbols.position
_theta_finite = atan2(_r, _x)
_theta_infinite = limit(_theta_finite, _x, S.Infinity)
assert _theta_infinite == 0

_magnetic_flux_density = _finite_wire_law.law.rhs.subs({
    _finite_wire_law.absolute_permeability: absolute_permeability,
    _finite_wire_law.current: current,
    _finite_wire_law.first_angle: _theta_infinite,
    _finite_wire_law.second_angle: _theta_infinite,
    _finite_wire_law.radial_distance: radial_distance,
})

assert expr_equals(_magnetic_flux_density, law.rhs)


@validate_input(
    absolute_permeability_=absolute_permeability,
    current_=current,
    radial_distance_=radial_distance,
)
@validate_output(magnetic_flux_density)
def calculate_magnetic_field(
    absolute_permeability_: Quantity,
    current_: Quantity,
    radial_distance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        absolute_permeability: absolute_permeability_,
        current: current_,
        radial_distance: radial_distance_,
    })
    return Quantity(result)


# UNIQUE_LAW_ID: 486
