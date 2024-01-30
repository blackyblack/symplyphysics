from sympy import Eq, solve, sqrt, S
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output, dimensionless, convert_to)


# Description
# Froude number characterizes the ratio between the force of inertia and
# the external force, in the field of which the motion occurs, acting on
# an elementary volume of liquid or gas.
# Law: Fr = u / sqrt(g * L), where
# u is velocity,
# g is acceleration due to gravity,
# L is a characteristic length.

velocity = Symbol("velocity", units.velocity)
characteristic_length = Symbol("characteristic_length", units.length)
froude_number = Symbol("froude_number", dimensionless)

law = Eq(froude_number, velocity /
         sqrt(units.acceleration_due_to_gravity * characteristic_length))


def print_law() -> str:
    return print_expression(law)


@validate_input(velocity_=velocity, characteristic_length_=characteristic_length)
@validate_output(froude_number)
def calculate_froude_number(velocity_: Quantity, characteristic_length_: Quantity) -> float:
    result_expr = solve(law, froude_number, dict=True)[0][froude_number]
    result_applied = result_expr.subs({
        velocity: velocity_,
        characteristic_length: characteristic_length_
    })
    result = Quantity(result_applied)
    return float(convert_to(result, S.One).evalf())
