"""
Weight in fluid via ratio of densities
======================================

The *Archimedean force* acting on a body fully submersed in a fluid is equal to the weight of the
fluid displaced by the body. It can be derived that the weight of the body submersed in the fluid
is proportional to its weight in vacuum and also depends on the ratio of the fluid density and body
density.

**Conditions:**

#. The body is completely *submersed* in the fluid.

#. The density of the body is no less than than the density of the fluid, otherwise the body would
   float to the surface and will be only partially submerged.

**Links:**

#. `Physics LibreTexts, derivable from here <https://phys.libretexts.org/Bookshelves/University_Physics/Physics_(Boundless)/10%3A_Fluids/10.3%3A_Archimedes_Principle>`__.
"""

from sympy import Eq, solve, pi, Idx, ask
from symplyphysics import (clone_as_symbol, symbols, Quantity, validate_input, validate_output,
    quantities, global_index)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.definitions import (
    density_from_mass_volume as _density_def,
    net_force_is_sum_of_individual_forces as _superposition_law,
)
from symplyphysics.laws.dynamics import (
    acceleration_is_force_over_mass as _newtons_law,
    buoyant_force_from_density_and_volume as _archimedes_law,
)
from symplyphysics.laws.geometry import (
    scalar_projection_is_vector_length_times_cosine_of_angle as _projection_law,)

weight_in_fluid = clone_as_symbol(symbols.force,
    display_symbol="W_fl",
    display_latex="W_\\text{fl}")
"""
Weight of the body submersed in the fluid. See :symbols:`force`.
"""

weight_in_vacuum = clone_as_symbol(symbols.force,
    display_symbol="W_vac",
    display_latex="W_\\text{vac}")
"""
Weight of the body in vacuum, i.e. its true weight. See :symbols:`force`.
"""

fluid_density = clone_as_symbol(
    symbols.density,
    display_symbol="rho_fl",
    display_latex="\\rho_\\text{fl}",
    positive=True,
)
"""
:symbols:`density` of the fluid.
"""

body_density = clone_as_symbol(
    symbols.density,
    display_symbol="rho_b",
    display_latex="\\rho_\\text{b}",
    positive=True,
)
"""
:symbols:`density` of the body.
"""

law = Eq(weight_in_fluid, weight_in_vacuum * (1 - (fluid_density / body_density)))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

_body_volume = clone_as_symbol(_density_def.volume, positive=True)

# 1. Find the Archimedes' force and project it onto the gravity force vector.

# Note that the buoyant force is co-directed with the z-axis, so this expression is also the
# projection of the buoyant force vector to the z-axis.
_archimedes_force_expr = _archimedes_law.law.rhs.subs({
    _archimedes_law.fluid_density: fluid_density,
    _archimedes_law.displaced_volume: _body_volume,
})

# 2. True weight of the body is equal to the gravity force on the body as a whole.

_body_mass_expr = solve(
    _density_def.definition,
    _density_def.mass,
)[0].subs({
    _density_def.density: body_density,
    _density_def.volume: _body_volume,
})

# Equal to the true weight of the body
_gravity_force_norm_expr = solve(
    _newtons_law.law,
    _newtons_law.force,
)[0].subs({
    _newtons_law.acceleration: quantities.acceleration_due_to_gravity,
    _newtons_law.mass: _body_mass_expr,
})

# The force of gravity is parallel to the vector of acceleration due to gravity, which in turn is
# parallel to the z-axis but points in the other direction, therefore the angle between the force
# of gravity and the unit vector parallel to the z-axis is `pi`
_gravity_force_proj_expr = _projection_law.law.rhs.subs({
    _projection_law.vector_length: _gravity_force_norm_expr,
    _projection_law.angle: pi,
})

# 3. Apparent weight is the true weight with the Archimedes' force substituted from it since the
#    gravity force and the buoyant force are antiparallel.

# This is the projection of the net force on the z-axis.
_net_force_proj_expr = _superposition_law.definition.rhs.subs(
    global_index,
    Idx("i", (1, 2)),
).doit().subs({
    _superposition_law.force[1]: _gravity_force_proj_expr,
    _superposition_law.force[2]: _archimedes_force_expr,
})

# The net force points opposite to direction of the z-axis.
assert ask(_net_force_proj_expr < 0, body_density > fluid_density)

# Therefore the apparent weight, being the absolute value of the net force on the body, can be
# found like this, note that `abs(x) = -1 * x` if `x < 0`.
_apparent_weight_expr = -1 * _net_force_proj_expr

# 4. Perform final replacements

_true_weight_eqn = Eq(weight_in_vacuum, _gravity_force_norm_expr)
_apparent_weight_eqn = Eq(weight_in_fluid, _apparent_weight_expr)

_apparent_weight_expr = solve(
    (_true_weight_eqn, _apparent_weight_eqn),
    (_body_volume, weight_in_fluid),
    dict=True,
)[0][weight_in_fluid]

assert expr_equals(_apparent_weight_expr, law.rhs)


@validate_input(weight_vacuum_=weight_in_vacuum,
    liquid_density_=fluid_density,
    body_density_=body_density)
@validate_output(weight_in_fluid)
def calculate_weight(weight_vacuum_: Quantity, liquid_density_: Quantity,
    body_density_: Quantity) -> Quantity:
    result_expr = solve(law, weight_in_fluid, dict=True)[0][weight_in_fluid]
    result_weight = result_expr.subs({
        weight_in_vacuum: weight_vacuum_,
        fluid_density: liquid_density_,
        body_density: body_density_,
    })

    return Quantity(result_weight)
