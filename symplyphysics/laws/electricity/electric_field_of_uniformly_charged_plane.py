r"""
Electric field of uniformly charged plane
=========================================

The electric field strength of a uniformly charged plane is proportional to its
charge density.

**Notation:**

#. :math:`\varepsilon_0` (:code:`epsilon_0`) is vacuum permittivity.

**Conditions:**

#. The plane is thin, i.e. its thickness approaches zero.
#. The plane is in vacuum.

**Links:**

#. `Electric field of a uniformly charged plane <https://farside.ph.utexas.edu/teaching/316/lectures/node27.html>`__
"""

from sympy import (Eq, solve, symbols as sympy_symbols)
from symplyphysics import (
    units,
    quantities,
    Quantity,
    SymbolNew,
    Vector,
    validate_input,
    validate_output,
    scale_vector,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity.vector import (
    electric_flux_of_uniform_electric_field as _flux_law,
)
from symplyphysics.laws.electricity import (
    electric_flux_through_closed_surface_via_total_charge as _gauss_law,
)
from symplyphysics.laws.quantities import (
    quantity_is_areal_density_times_area as _areal_qty_law,
)

electric_field_strength = SymbolNew("E", units.voltage / units.length)
"""
Value of the electric field.
"""

surface_charge_density = SymbolNew("sigma", units.charge / units.area, display_latex="\\sigma")
"""
Surface charge density of the plane.
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

_electric_field_left = Vector([-1 * electric_field_strength, 0, 0])
_electric_field_right = Vector([electric_field_strength, 0, 0])

_cross_sectional_area = sympy_symbols("cross_sectional_area")

_area_left = scale_vector(_cross_sectional_area, Vector([-1, 0, 0]))
_area_right = scale_vector(_cross_sectional_area, Vector([1, 0, 0]))

_electric_flux_left = _flux_law.electric_flux_law(
    electric_field_=_electric_field_left,
    area_=_area_left,
)

_electric_flux_right = _flux_law.electric_flux_law(
    electric_field_=_electric_field_right,
    area_=_area_right,
)

# The electric field is orthogonal to the normal vector of the cylinder's side at all points.
# This, the flux there would be zero.
_electric_flux_side = 0

# The whole integration area is composed of the two cross-sections and the cylinder side.
_total_electric_flux = _electric_flux_left + _electric_flux_right + _electric_flux_side

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
