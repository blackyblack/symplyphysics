from sympy import Eq, solve, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics.springs import spring_reaction_is_proportional_to_deformation as hookes_law

# Description
## If two springs are side-to-side to one another, i.e. connected in parallel, the total
## stiffness of the system of springs is equal to the sum of the stiffnesses of each spring.

# Law: k = k1 + k2
## k - total stiffness
## k1 - stiffness of first spring
## k2 - stiffness of second spring

# Condition
## - Springs must be Hookean, or linear-response, i.e. obey the Hooke's law.

total_stiffness = Symbol("total_stiffness", units.force / units.length)
first_stiffness = Symbol("first_stiffness", units.force / units.length)
second_stiffness = Symbol("second_stiffness", units.force / units.length)

law = Eq(total_stiffness, first_stiffness + second_stiffness)

# Derive law from Hooke's law
## In the case of serial connection, the forces acting on both springs add up while
## the springs deform equally.

_deformation = SymSymbol("deformation")

_force_expr = solve(
    hookes_law.law, hookes_law.spring_reaction
)[0].subs(
    hookes_law.deformation, _deformation
)

_first_force = _force_expr.subs(hookes_law.stiffness, first_stiffness)
_second_force = _force_expr.subs(hookes_law.stiffness, second_stiffness)

_total_force_hooke = _force_expr.subs(hookes_law.stiffness, total_stiffness)
_total_force_added = _first_force + _second_force

_total_stiffness_derived = solve(
    Eq(_total_force_hooke, _total_force_added),
    total_stiffness,
)[0]

assert expr_equals(_total_stiffness_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    first_stiffness_=first_stiffness,
    second_stiffness_=second_stiffness,
)
@validate_output(total_stiffness)
def calculate_total_stiffness(
    first_stiffness_: Quantity,
    second_stiffness_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        first_stiffness: first_stiffness_,
        second_stiffness: second_stiffness_,
    })
    return Quantity(result)
