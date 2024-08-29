from sympy import (Eq, solve, pi)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## A magnetic field intensity is a vector field that describes the magnetic influence on moving electric
## charges, electric currents and magnetic materials. In the case of an infinite wire, the magnetic
## field intensity depends only on the current and the distance from the wire.

## Law is: H = I / (2 * pi * r), where
## H - magnetic field intensity,
## I - current,
## r - distance from the wire.

magnetic_intensity = Symbol("magnetic_intensity", units.current / units.length)

current = Symbol("current", units.current)
distance = Symbol("distance", units.length)

law = Eq(magnetic_intensity, current / (2 * pi * distance))


@validate_input(current_=current, distance_=distance)
@validate_output(magnetic_intensity)
def calculate_magnetic_intensity(current_: Quantity, distance_: Quantity) -> Quantity:
    result_expr = solve(law, magnetic_intensity, dict=True)[0][magnetic_intensity]
    result_expr = result_expr.subs({
        current: current_,
        distance: distance_,
    })
    return Quantity(result_expr)
