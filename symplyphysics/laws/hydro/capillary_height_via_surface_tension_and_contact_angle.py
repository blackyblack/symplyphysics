"""
Capillary height via surface tension and contact angle
=======================================================

The **Jurin's law** determines the height to which the liquid rises in capillaries. It
states that the maximum height of a liquid in a capillary tube is inversely proportional
to the tube's diameter.

**Notation:**

#. :quantity_notation:`acceleration_due_to_gravity`.

**Conditions:**

#. The surface of the meniscus is spherical.
#. Height :math:`h` of the raised (lowered) liquid is much larger than the radius
   :math:`r` of the capillary.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Jurin%27s_law>`__.

..
    TODO: rename file to use descriptive name
"""

from sympy import Eq, solve, cos, dsolve
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities,
    clone_as_symbol)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.thermodynamics import (
    laplace_pressure_of_spherical_shapes as _young_laplace_law,)
from symplyphysics.laws.hydro import (
    hydrostatic_pressure_via_density_and_height as _hydrostatic_pressure_law,
    inner_pressure_is_constant as _bernoulli_principle,
    inner_pressure_is_sum_of_pressures as _inner_pressure_def,
    laplace_pressure_is_pressure_difference as _laplace_pressure_def,
)

height = symbols.height
"""
:symbols:`height` of the liquid column.
"""

surface_tension = symbols.surface_tension
"""
:symbols:`surface_tension` of the liquid.
"""

angle = symbols.angle
"""
Contact :symbols:`angle` between of the liquid and the tube wall.
"""

density = symbols.density
"""
:symbols:`density` of the liquid.
"""

radius = symbols.radius
"""
:symbols:`radius` of the capillary.
"""

law = Eq(
    height,
    2 * surface_tension * cos(angle) / (density * radius * quantities.acceleration_due_to_gravity))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

_laplace_pressure = _young_laplace_law.laplace_pressure
_inner_pressure = _bernoulli_principle.inner_pressure
_time = _bernoulli_principle.time
_height = _hydrostatic_pressure_law.height
_atmospheric_pressure = clone_as_symbol(symbols.pressure, subscript="atm")

# The meniscus has a spherical shape (since a sphere is a shape with the smallest surface area for
# a given volume). Also see [this figure](https://en.wikipedia.org/wiki/Young%E2%80%93Laplace_equation#/media/File:Spherical_meniscus.PNG)
_radius_of_curvature = radius / cos(angle)

_young_laplace_eqn = _young_laplace_law.law.subs({
    _young_laplace_law.surface_tension: surface_tension,
    _young_laplace_law.radius_of_curvature: _radius_of_curvature,
})

_hydrostatic_pressure_expr = _hydrostatic_pressure_law.law.rhs.subs({
    _hydrostatic_pressure_law.density: density,
})

_bernoulli_dsolved_eqn = dsolve(
    _bernoulli_principle.law,
    _inner_pressure(_time),
    ics={_inner_pressure(0): _inner_pressure(0)},
)

_inner_pressure_expr = _inner_pressure_def.law.rhs.subs({
    _inner_pressure_def.dynamic_pressure: 0,  # the fluid is static
    _inner_pressure_def.hydrostatic_pressure: _hydrostatic_pressure_expr,
    _inner_pressure_def.inner_pressure: _inner_pressure(_time),
})

# In the following computations, both "outside inner pressure" and "inside inner pressure" are
# measured at equal heights.

# a. Concave meniscus: the interface curves towards the inside of the tube. The inside pressure is
#    the inner pressure under the meniscus, and the outside pressure is the atmospheric pressure.

# Outside the tube at the fluid's surface
_outside_inner_pressure_expr = _inner_pressure_expr.subs({
    _time: 0,
    _height: 0,
    _inner_pressure_def.static_pressure: _atmospheric_pressure,
})

# Inside the tube at the bottom of the fluid column, which is the same height as the fluid's
# surface outside the tube.
_inside_inner_pressure_expr = _inner_pressure_expr.subs({
    _inner_pressure_def.static_pressure: _laplace_pressure_def.inside_pressure,
    _height: height,
})

_laplace_pressure_eqn = _laplace_pressure_def.law.subs({
    _laplace_pressure_def.laplace_pressure: _laplace_pressure,
    _laplace_pressure_def.outside_pressure: _atmospheric_pressure,
})

_bernoulli_eqn = _bernoulli_dsolved_eqn.subs({
    _inner_pressure(0): _outside_inner_pressure_expr,
    _inner_pressure(_time): _inside_inner_pressure_expr,
})

_height_derived = solve(
    (_young_laplace_eqn, _bernoulli_eqn, _laplace_pressure_eqn),
    (height, _laplace_pressure, _laplace_pressure_def.inside_pressure),
    dict=True,
)[0][height]

assert expr_equals(_height_derived, law.rhs)

# b. Convex meniscus: the interface curves towards the outside of the tube (i.e. the atmosphere).
#    The inside pressure is the atmospheric pressure and the outside pressure is the pressure under
#    the meniscus.

# Outside of the tube at a `height` (below the surface) of the lowered fluid in the tube.
_outside_inner_pressure_expr = _inner_pressure_expr.subs({
    _time: 0,
    _height: height,
    _inner_pressure_def.static_pressure: _atmospheric_pressure,
})

# Inside the tube at the surface of the meniscus, which is lowered at a certain `height` relative
# to the fluid's surface outside of the meniscus.
_inside_inner_pressure_expr = _inner_pressure_expr.subs({
    _inner_pressure_def.static_pressure: _laplace_pressure_def.outside_pressure,
    _height: 0,
})

_laplace_pressure_eqn = _laplace_pressure_def.law.subs({
    _laplace_pressure_def.laplace_pressure: _laplace_pressure,
    _laplace_pressure_def.inside_pressure: _atmospheric_pressure,
})

_bernoulli_eqn = _bernoulli_dsolved_eqn.subs({
    _inner_pressure(0): _outside_inner_pressure_expr,
    _inner_pressure(_time): _inside_inner_pressure_expr,
})

_height_derived = solve(
    (_young_laplace_eqn, _bernoulli_eqn, _laplace_pressure_eqn),
    (height, _laplace_pressure, _laplace_pressure_def.outside_pressure),
    dict=True,
)[0][height]

assert expr_equals(_height_derived, law.rhs)


@validate_input(surface_tension_coefficient_=surface_tension,
    angle_=angle,
    density_of_liquid_=density,
    radius_=radius)
@validate_output(height)
def calculate_height(surface_tension_coefficient_: Quantity, angle_: Quantity,
    density_of_liquid_: Quantity, radius_: Quantity) -> Quantity:
    result_expr = solve(law, height, dict=True)[0][height]
    result_height = result_expr.subs({
        surface_tension: surface_tension_coefficient_,
        angle: angle_,
        density: density_of_liquid_,
        radius: radius_,
    })
    result_height_quantity = Quantity(result_height)
    if result_height_quantity.scale_factor < radius_.scale_factor:
        raise ValueError(
            f"The height must be greater than the radius. Currently {result_height_quantity.scale_factor} < {radius_.scale_factor}"
        )
    return result_height_quantity
