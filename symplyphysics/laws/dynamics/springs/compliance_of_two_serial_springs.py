"""
Compliance of two serial springs
================================

If two springs are connected end-to-end to one another, i.e. connected in series, the
total compliance of the system of springs is the sum of the compliances of each spring.

**Conditions:**

#. Springs must be Hookean, or linear-response, i.e. obey the Hooke's law.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Series_and_parallel_springs#Formulas>`__.
"""

from sympy import Eq, solve, Symbol as SymSymbol
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import compliance_is_inverse_stiffness as compliance_def
from symplyphysics.laws.dynamics import reaction_force_from_action_force as newtons_third_law
from symplyphysics.laws.dynamics.springs import spring_reaction_is_proportional_to_deformation as hookes_law

total_compliance = symbols.compliance
"""
Total :symbols:`compliance`.
"""

first_compliance = clone_as_symbol(symbols.compliance, subscript="1")
"""
:symbols:`compliance` of the first spring.
"""

second_compliance = clone_as_symbol(symbols.compliance, subscript="2")
"""
:symbols:`compliance` of the second spring.
"""

law = Eq(total_compliance, first_compliance + second_compliance)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from Hooke's law
## In the case of serial connection, the forces acting on both springs are equal whereas
## their deformations add up.

# When external force is applied to the first spring, it produces a reaction force on both sides of
# that spring. That in turn produces a reaction force on both sides of the other spring.
_external_force = SymSymbol("_external_force")
_first_spring_reaction = abs(
    newtons_third_law.law.rhs.subs(newtons_third_law.action_force, _external_force))
_second_spring_reaction = abs(
    newtons_third_law.law.rhs.subs(newtons_third_law.action_force, _first_spring_reaction))
_total_reaction = abs(
    newtons_third_law.law.rhs.subs(newtons_third_law.action_force, _second_spring_reaction))

_stiffness_expr = solve(compliance_def.definition, compliance_def.stiffness)[0]
_first_stiffness = _stiffness_expr.subs(compliance_def.compliance, first_compliance)
_second_stiffness = _stiffness_expr.subs(compliance_def.compliance, second_compliance)
_total_stiffness = _stiffness_expr.subs(compliance_def.compliance, total_compliance)

_deformation_expr = solve(hookes_law.law, hookes_law.deformation)[0]

_first_deformation = _deformation_expr.subs({
    hookes_law.spring_reaction: _first_spring_reaction,
    hookes_law.stiffness: _first_stiffness,
})
_second_deformation = _deformation_expr.subs({
    hookes_law.spring_reaction: _second_spring_reaction,
    hookes_law.stiffness: _second_stiffness,
})

_total_deformation_hooke = _deformation_expr.subs({
    hookes_law.spring_reaction: _total_reaction,
    hookes_law.stiffness: _total_stiffness,
})

_total_deformation_added = _first_deformation + _second_deformation

_total_compliance_derived = solve(
    Eq(_total_deformation_hooke, _total_deformation_added),
    total_compliance,
)[0]

assert expr_equals(_total_compliance_derived, law.rhs)


@validate_input(
    first_compliance_=first_compliance,
    second_compliance_=second_compliance,
)
@validate_output(total_compliance)
def calculate_compliance(
    first_compliance_: Quantity,
    second_compliance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        first_compliance: first_compliance_,
        second_compliance: second_compliance_,
    })
    return Quantity(result)
