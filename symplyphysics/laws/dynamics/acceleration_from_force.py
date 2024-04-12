from sympy import (Eq, solve, sympify)
from symplyphysics import (Vector, Quantity, print_expression, validate_input, validate_output,
    symbols)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics.vector import acceleration_from_force as acceleration_law_vector

# Description
## Newton's second law: a = F / m

acceleration = symbols.kinematic.acceleration
force = symbols.dynamics.force
mass = symbols.basic.mass

law = Eq(acceleration, force / mass)

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


@validate_input(mass_=mass, acceleration_=acceleration)
@validate_output(force)
def calculate_force(mass_: Quantity, acceleration_: Quantity) -> Quantity:
    result_force_expr = solve(law, force, dict=True)[0][force]
    result_expr = result_force_expr.subs({
        mass: mass_,
        acceleration: acceleration_
    })
    return Quantity(result_expr)
