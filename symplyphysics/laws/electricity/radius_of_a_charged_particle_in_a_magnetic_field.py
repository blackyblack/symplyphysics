from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Let an arbitrary particle move in a magnetic field around a circle. Then radius of its motion depends on
## magnetic induction and on velocity, charge and mass of the particle.

## Law is: r = m * v / (q * B), where
## r - radius,
## m - mass,
## v - velocity,
## q - charge,
## B - induction.

radius = Symbol("radius", units.length)

mass = Symbol("mass", units.mass)
velocity = Symbol("velocity", units.velocity)
charge = Symbol("charge", units.charge)
induction = Symbol("induction", units.magnetic_density)

law = Eq(radius, mass * velocity / (charge * induction))


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=mass, velocity_=velocity, induction_=induction, charge_=charge)
@validate_output(radius)
def calculate_radius(mass_: Quantity, velocity_: Quantity, induction_: Quantity, charge_: Quantity) -> Quantity:
    result_expr = solve(law, radius, dict=True)[0][radius]
    result_expr = result_expr.subs({
        mass: mass_,
        velocity: velocity_,
        induction: induction_,
        charge: charge_,
    })
    return Quantity(result_expr)
