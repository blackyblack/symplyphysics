from sympy import (Eq, solve, sympify)
from symplyphysics import (Vector, Quantity, print_expression, validate_input, validate_output,
    symbols)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics.vector import acceleration_from_force as acceleration_law_vector

# Description
## Newton's second law: a = F / m

law = Eq(symbols.kinematic.acceleration, symbols.dynamics.force / symbols.basic.mass)

# Derive the same law from vector form

# Scalar law is equivalent to using one-dimensional vectors
force_vector = Vector([symbols.dynamics.force])
acceleration_vector = acceleration_law_vector.acceleration_law(force_vector)
assert len(acceleration_vector.components) == 1
acceleration_with_mass = sympify(acceleration_vector.components[0]).subs(
    acceleration_law_vector.symbols.basic.mass, symbols.basic.mass)
assert expr_equals(acceleration_with_mass, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=symbols.basic.mass, acceleration_=symbols.kinematic.acceleration)
@validate_output(symbols.dynamics.force)
def calculate_force(mass_: Quantity, acceleration_: Quantity) -> Quantity:
    result_force_expr = solve(law, symbols.dynamics.force, dict=True)[0][symbols.dynamics.force]
    result_expr = result_force_expr.subs({
        symbols.basic.mass: mass_,
        symbols.kinematic.acceleration: acceleration_
    })
    return Quantity(result_expr)
