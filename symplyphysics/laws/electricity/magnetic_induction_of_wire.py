from sympy import (Eq, solve, pi)
from sympy.physics.units import magnetic_constant
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, dimensionless)

# Description
## Let there be an infinite thin conductor. Then the magnetic field created by the current in
## the conductor will depend only on the current, the material and the distance from the conductor.

## Law is: B = mu0 * mu * I / (2 * pi * r), where
## B - magnetic field induction,
## mu0 - magnetic_constant,
## mu - relative permeability of medium,
## I - current,
## r - distance from conductor.

# Conditions:
## - The wire is infinitely thin.

induction = Symbol("induction", units.magnetic_density)

relative_permeability = Symbol("relative_permeability", dimensionless)
current = Symbol("current", units.current)
distance = Symbol("distance", units.length)

law = Eq(induction, magnetic_constant * relative_permeability * current / (2 * pi * distance))


@validate_input(relative_permeability_=relative_permeability, current_=current, distance_=distance)
@validate_output(induction)
def calculate_induction(relative_permeability_: float, current_: Quantity,
    distance_: Quantity) -> Quantity:
    result_expr = solve(law, induction, dict=True)[0][induction]
    result_expr = result_expr.subs({
        relative_permeability: relative_permeability_,
        current: current_,
        distance: distance_
    })
    return Quantity(result_expr)
