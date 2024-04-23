from sympy import (Eq, solve)
from symplyphysics import (clone_symbol, symbols, Quantity, Symbol, print_expression, dimensionless,
    validate_input, validate_output)

# Description
## Friction force is tangential interaction between two objects, which impedes there relative movement.
## It is proportional to pressure of one object to another.

## Law: Ffriction = mu * N
## Where:
## Ffriction is friction force,
## mu is friction factor for this pair of objects,
## N is normal reaction from one object to another.

friction_force = clone_symbol(symbols.dynamics.force, "friction_force")
friction_factor = Symbol("friction_factor", dimensionless)
normal_reaction = clone_symbol(symbols.dynamics.force, "normal_reaction")

law = Eq(friction_force, friction_factor * normal_reaction)


def print_law() -> str:
    return print_expression(law)


@validate_input(friction_factor_=friction_factor, normal_reaction_=normal_reaction)
@validate_output(friction_force)
def calculate_friction_force(friction_factor_: float, normal_reaction_: Quantity) -> Quantity:
    result_expr = solve(law, friction_force, dict=True)[0][friction_force]
    friction_force_applied = result_expr.subs({
        friction_factor: friction_factor_,
        normal_reaction: normal_reaction_
    })
    return Quantity(friction_force_applied)
