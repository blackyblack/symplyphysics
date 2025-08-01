"""
Electric field of uniformly charged plane
=========================================

The electric field strength of a uniformly charged plane is proportional to its
charge density.

**Notation:**

#. :quantity_notation:`vacuum_permittivity`.

**Conditions:**

#. The plane is thin, i.e. its thickness approaches zero.
#. The plane is in vacuum.

**Links:**

#. `Electric field of a uniformly charged plane <https://farside.ph.utexas.edu/teaching/316/lectures/node27.html>`__
"""

from sympy import Eq, solve, symbols as sympy_symbols
from symplyphysics import quantities, Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity.vector import (
    electric_flux_of_uniform_electric_field as _flux_law,)
from symplyphysics.laws.electricity import (
    electric_flux_through_closed_surface_via_total_charge as _gauss_law,)
from symplyphysics.laws.quantities import (
    quantity_is_areal_density_times_area as _areal_qty_law,)

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, CoordinateVector

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

surface_charge_density = symbols.surface_charge_density
"""
:symbols:`surface_charge_density`.
"""

law = Eq(electric_field_strength, surface_charge_density / (2 * quantities.vacuum_permittivity))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from Gauss's law

# Imagine a cylinder as a gaussian surface whose axis is parallel to the electric field due
# to the charged plane, which we can set to be parallel to the x-axis without loss of
# generality. It is positioned in such a way that the plane cuts it in two equal halves.
# Let `a` be half of the axial dimension of the cylinder and `S` its cross-sectional area.
# Then for the projections of electric field `E` we would have `E(a) = -1 * E(-1 * a)` since
# the direction of the electric field in one half-space is opposite to that in the other half-space.

_cross_sectional_area = sympy_symbols("cross_sectional_area")

_electric_field_left = CoordinateVector([-1 * electric_field_strength, 0, 0], CARTESIAN)
_electric_field_right = CoordinateVector([electric_field_strength, 0, 0], CARTESIAN)

_e_x = CoordinateVector([1, 0, 0], CARTESIAN)

# The area pseudovector has the magnitude of the area of the surface and the direction of the unit
# normal to the surface. It is needed for the proper calculation of the electric flux. Also see
# `area` in `./vector/electric_flux_of_uniform_electric_field`
_vector_area_left = -1 * _cross_sectional_area * _e_x
_vector_area_right = _cross_sectional_area * _e_x

_electric_flux_left = _flux_law.law.rhs.subs({
    _flux_law.electric_field: _electric_field_left,
    _flux_law.area: _vector_area_left,
}).doit()

_electric_flux_right = _flux_law.law.rhs.subs({
    _flux_law.electric_field: _electric_field_right,
    _flux_law.area: _vector_area_right,
}).doit()

# The electric field is orthogonal to the normal vector of the cylinder's side at all points.
# This, the flux there would be zero.
ELECTRIC_FLUX_SIDE = 0

# The whole integration area is composed of the two cross-sections and the cylinder side.
_total_electric_flux = _electric_flux_left + _electric_flux_right + ELECTRIC_FLUX_SIDE

# The total charge of the cylinder is contained in the part of the charged plane that is
# contained within the cylinder.
_total_charge = _areal_qty_law.law.rhs.subs({
    _areal_qty_law.areal_density: surface_charge_density,
    _areal_qty_law.area: _cross_sectional_area,
})

_gauss_eqn = _gauss_law.law.subs({
    _gauss_law.total_electric_flux: _total_electric_flux,
    _gauss_law.total_charge: _total_charge,
})

_electric_field_expr = solve(_gauss_eqn, electric_field_strength)[0]

assert expr_equals(_electric_field_expr, law.rhs)


@validate_input(surface_charge_density_=surface_charge_density)
@validate_output(electric_field_strength)
def calculate_electric_intensity(surface_charge_density_: Quantity) -> Quantity:
    result_expr = law.rhs.subs({
        surface_charge_density: surface_charge_density_,
    })
    return Quantity(result_expr)
