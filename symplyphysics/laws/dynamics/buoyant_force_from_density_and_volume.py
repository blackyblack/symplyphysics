"""
Buoyant force from density and volume
=====================================

Any object, totally or partially immersed in a fluid (i.e. liquid or gas), is buoyed up by a force equal to the
weight of the fluid displaced by the object. Also known as the Archimedes principle. The *buoyant force*
vector is directed opposite to the gravity vector.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Archimedes%27_principle#Formula>`__.
"""

from sympy import Eq, solve, Idx
from symplyphysics import (clone_as_symbol, symbols, quantities, Quantity, validate_input,
    validate_output, global_index)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.core.experimental.coordinate_systems import CoordinateVector, CARTESIAN
from symplyphysics.core.experimental.vectors import VectorSymbol, VectorDot

from symplyphysics.definitions.vector import superposition_of_forces_is_sum as _superposition_law
from symplyphysics.laws.hydro import hydrostatic_pressure_via_density_and_height as _pressure_depth_law

buoyant_force = clone_as_symbol(symbols.force, display_symbol="F_A", display_latex="F_\\text{A}")
"""
The buoyant (Archimedes) :symbols:`force`.
"""

fluid_density = symbols.density
"""
The :symbols:`density` of the fluid.
"""

displaced_volume = symbols.volume
"""
The :symbols:`volume` of the displaced fluid. Equivalently, the volume of the part of the body
immersed in the fluid.
"""

law = Eq(buoyant_force,
    -1 * fluid_density * quantities.acceleration_due_to_gravity * displaced_volume)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law for a parallelepiped-shaped object whose sides are aligned with the axes of the
# Cartesian coordinate system.

_z_length = clone_as_symbol(symbols.length, subscript="z")
_xy_area = clone_as_symbol(symbols.area, subscript="xy")

# Depth is the difference in the `z` coordinates between the given point and the top side
_depth = _pressure_depth_law.height

_top_depth = clone_as_symbol(_depth, subscript="0")
_bottom_depth = _top_depth + _z_length

# a.k.a. gauge pressure; it is a function of `_depth`
_hydrostatic_pressure = _pressure_depth_law.law.rhs.subs(_pressure_depth_law.density, fluid_density)

_unit_normal = VectorSymbol("n")
_area = symbols.area

# Note that we ignore the external (e.g. atmospheric) pressure because it is considered constant
# and it therefore does not contribute to the buoyant force.

# TODO: use law here
_force_expr = -1 * _hydrostatic_pressure * _unit_normal * _area

_k = CoordinateVector([0, 0, 1], CARTESIAN)

_top_force_expr = _force_expr.subs({
    _depth: _top_depth,
    _unit_normal: _k,
    _area: _xy_area,
})

_bottom_force_expr = _force_expr.subs({
    _depth: _bottom_depth,
    _unit_normal: -1 * _k,
    _area: _xy_area,
})

# NOTE: the forces exerted on the sides cancel each other out. For instance, the 2 sides parallel
# to the xOz plane experience forces of the same magnitude (since the force only depends on the
# depth) but opposite direction; the same applies to the sides parallel to the yOz plane. But to
# show this here, vector integration should be implemented.

_net_force_expr = _superposition_law.law.rhs.subs(
    global_index,
    Idx("i", (1, 2)),
).doit().subs({
    _superposition_law.force[1]: _top_force_expr,
    _superposition_law.force[2]: _bottom_force_expr,
})
_net_force_expr = CoordinateVector.from_expr(_net_force_expr)

# The axis along which depth is measured in this law is antiparallel to the z-axis
_net_force_projection = VectorDot(_net_force_expr, -1 * _k).doit()

# Replace `l_z * A_xy` with `V`
_net_force_projection = solve(
    [
    Eq(buoyant_force, _net_force_projection),
    Eq(displaced_volume, _z_length * _xy_area),
    ],
    (buoyant_force, _z_length),
    dict=True,
)[0][buoyant_force]

assert expr_equals(law.rhs, _net_force_projection)

# For any other shape of the object, we can approximate it with parallelepipeds


@validate_input(fluid_density_=fluid_density, displaced_volume_=displaced_volume)
@validate_output(buoyant_force)
def calculate_force_buoyant(fluid_density_: Quantity, displaced_volume_: Quantity) -> Quantity:
    result_force_expr = solve(law, buoyant_force, dict=True)[0][buoyant_force]
    result_expr = result_force_expr.subs({
        fluid_density: fluid_density_,
        displaced_volume: displaced_volume_
    })
    return Quantity(abs(result_expr))
