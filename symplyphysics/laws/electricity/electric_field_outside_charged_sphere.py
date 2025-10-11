"""
Electric field outside charged sphere
=====================================

The electric field outside of a charged sphere behaves as is the sphere were a point charge,
i.e. its magnitude is inversely proportional to the square of the distance to the center
of the sphere. However, due to the Gauss's law, on the inside the electric field is exactly zero.

**Notation:**

#. :quantity_notation:`vacuum_permittivity`.

**Conditions:**

#. The sphere is thin, i.e. its thickness approaches zero.

#. The medium is vacuum.

**Links:**

#. `Wikipedia, "Spherical volume" <https://en.wikipedia.org/wiki/Electric_field#Common_formul%C3%A6>`__.
"""

from sympy import (Eq, solve, pi, integrate)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities,
    clone_as_symbol)
from symplyphysics.definitions.vector import vector_area_is_unit_normal_times_scalar_area as _vector_area_def
from symplyphysics.laws.electricity import electric_flux_through_closed_surface_via_total_charge as _gauss_law
from symplyphysics.laws.electricity.vector import electric_flux_of_uniform_electric_field as _flux_def

from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.symbols import BasicSymbol
from symplyphysics.core.experimental.coordinate_systems import SPHERICAL, CoordinateVector, AppliedPoint
from symplyphysics.core.experimental.coordinate_systems.surface import Surface

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

charge = symbols.charge
"""
Total :symbols:`charge` of the sphere.
"""

distance = clone_as_symbol(symbols.euclidean_distance, positive=True)
"""
:symbols:`euclidean_distance` to the center of the sphere.
"""

law = Eq(electric_field_strength, charge / (4 * pi * quantities.vacuum_permittivity * distance**2))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from Gauss's law, see :ref:`Electric flux through closed surface via total charge`

# Imagine a sphere `S'` concentric with the given sphere `S`. Imagine a small area `A` around a
# point `P` on `S'` and let us calculate the electric field flux through that area. The area
# vector `A` also has only a radial component (it points outside the sphere in the direction of
# the `r`-axis):

_p = BasicSymbol("P")
_e_r = CoordinateVector([1, 0, 0], SPHERICAL, _p)  # unit vector along the `r`-axis

# The `E`-field remains the same under any rotation whose axis passes through the origin (i.e.
# there is no preferred direction), therefore its components can only depend on `r` (the distance
# to origin). The rotational symmetry also means that the electric field only has the component
# along the `r` axis. So at any point of the sphere `S'`, the electric field has the same magnitude
# and is parallel to the radial vector:
_electric_field = electric_field_strength * _e_r

_, _theta, _phi = SPHERICAL.base_scalars
_sphere = Surface((_theta, _phi), AppliedPoint([distance, _theta, _phi], SPHERICAL))
# We do not multiply this by `d(theta) * d(phi)` because these differentials will be set to `1`
# when we integrate over the infinitesimal flux later on.
_area_change = _sphere.scalar_area_multiple

_vector_area_change = _vector_area_def.law.rhs.subs({
    _vector_area_def.scalar_area: _area_change,
    _vector_area_def.unit_normal: _e_r,
})

# The area is small enough for the electric field to be constant throughout it.
_flux_change = _flux_def.law.rhs.subs({
    _flux_def.electric_field: _electric_field,
    _flux_def.area: _vector_area_change,
}).doit()

# Now we can integrate the infinitesimal flux over the whole sphere:
_total_flux = integrate(_flux_change, (_theta, 0, pi), (_phi, 0, 2 * pi))

# Note that the radius of the sphere of integration must be greater than the radius of the charged
# sphere, otherwise `charge` would be 0.
_gauss_eqn = _gauss_law.law.subs({
    _gauss_law.total_electric_flux: _total_flux,
    _gauss_law.total_charge: charge,
})

_electric_field_strength_expr = solve(_gauss_eqn, electric_field_strength)[0]

assert expr_equals(_electric_field_strength_expr, law.rhs)


@validate_input(charge_=charge, distance_=distance)
@validate_output(electric_field_strength)
def calculate_electric_intensity(charge_: Quantity, distance_: Quantity) -> Quantity:
    result_expr = solve(law, electric_field_strength, dict=True)[0][electric_field_strength]
    result_expr = result_expr.subs({
        charge: charge_,
        distance: distance_,
    })
    return Quantity(result_expr)
