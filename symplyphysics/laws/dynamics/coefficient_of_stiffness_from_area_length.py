"""
Coefficient of stiffness from area and length
=============================================

The *stiffness* coefficient depends on the Young's modulus, cross-sectional area and length.
Young's modulus is a tabular value that is different for each material.

**Links:**

#. `Wikipedia, first formula <https://en.wikipedia.org/wiki/Stiffness#Relationship_to_elasticity>`__.
"""

from sympy import (Eq, solve, Q, Symbol as SymSymbol)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.dynamics import pressure_from_force_and_area as _pressure_law
from symplyphysics.laws.dynamics.springs import (
    spring_reaction_is_proportional_to_deformation as _hookes_law,)
from symplyphysics.laws.dynamics.deformation import (
    engineering_normal_strain_is_total_deformation_over_initial_dimension as _strain_def,
    tensile_stress_is_youngs_modulus_times_strain as _tensile_stress_law,
)

stiffness = symbols.stiffness
"""
The :symbols:`stiffness` coefficient of the material.
"""

young_modulus = symbols.young_modulus
"""
The :symbols:`young_modulus` of the material.
"""

area = symbols.area
"""
The cross-sectional :symbols:`area` of the object.
"""

length = symbols.length
"""
The :symbols:`length` of the object.
"""

law = Eq(stiffness, young_modulus * area / length)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law


def _replace(expr, eqn, *, old):
    """Replaces `old` within `expr` according to equation `eqn`."""

    sym = SymSymbol("x")
    solved = solve((Eq(sym, expr), eqn), (sym, old), dict=True)[-1][sym]
    return solved


# We're interested in the magnitude of the force, not its projection along the displacement vector
_abs_force_expr = abs(solve(_hookes_law.law, _hookes_law.spring_reaction)[0]).refine(
    Q.positive(_hookes_law.stiffness) & Q.positive(_hookes_law.deformation))

_hookes_eqn = Eq(_hookes_law.spring_reaction, _abs_force_expr)

_stiffness_expr = solve(_hookes_eqn, _hookes_law.stiffness)[0]

_pressure_eqn = _pressure_law.law.subs({
    _pressure_law.force: _hookes_law.spring_reaction,
    _pressure_law.area: area
})

_stiffness_expr = _replace(_stiffness_expr, _pressure_eqn, old=_hookes_law.spring_reaction)

_stress_eqn = _tensile_stress_law.law.subs({
    _tensile_stress_law.stress: _pressure_law.pressure,
    _tensile_stress_law.young_modulus: young_modulus,
})

_stiffness_expr = _replace(_stiffness_expr, _stress_eqn, old=_pressure_law.pressure)

_strain_eqn = _strain_def.law.subs({
    _strain_def.total_deformation: _hookes_law.deformation,
    _strain_def.initial_dimension: length,
})

_stiffness_expr = _replace(
    _stiffness_expr,
    _strain_eqn,
    old=_tensile_stress_law.engineering_normal_strain,
)

assert expr_equals(law.rhs, _stiffness_expr)


@validate_input(module_of_young_=young_modulus, area_=area, length_=length)
@validate_output(stiffness)
def calculate_coefficient_of_stiffness(module_of_young_: Quantity, area_: Quantity,
    length_: Quantity) -> Quantity:
    result_expr = solve(law, stiffness, dict=True)[0][stiffness]
    result_expr = result_expr.subs({young_modulus: module_of_young_, area: area_, length: length_})
    return Quantity(result_expr)
