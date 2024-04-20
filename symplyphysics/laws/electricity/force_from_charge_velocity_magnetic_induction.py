from sympy import (Eq, solve, sin)
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol, print_expression, validate_input,
    validate_output, angle_type)

# Description
## Lorentz force is force acting on a charge moving at speed from magnetic field.
## Lorentz force depends on the magnitude of charge, its velocity, magnitude of magnetic induction
## and angle between the magnetic induction and the charge velocity.
## The force is directed perpendicular to the plane in which the velocity and magnetic induction vectors are located.

## Law is: F = q * v * B * sin(a), where
## F - force,
## q - charge,
## v - velocity of charge,
## B - induction,
## a - angle between induction and velocity.

lorentz_force = clone_symbol(symbols.dynamics.force, "lorentz_force")

charge = Symbol("charge", units.charge)
velocity = Symbol("velocity", units.velocity)
induction = Symbol("induction", units.magnetic_density)
angle = Symbol("angle", angle_type)

law = Eq(lorentz_force, charge * velocity * induction * sin(angle))


def print_law() -> str:
    return print_expression(law)


@validate_input(charge_=charge, velocity_=velocity, angle_=angle, induction_=induction)
@validate_output(lorentz_force)
def calculate_force(charge_: Quantity, velocity_: Quantity, angle_: float | Quantity,
    induction_: Quantity) -> Quantity:
    result_expr = solve(law, lorentz_force, dict=True)[0][lorentz_force]
    result_expr = result_expr.subs({
        charge: charge_,
        velocity: velocity_,
        angle: angle_,
        induction: induction_,
    })
    return Quantity(result_expr)
