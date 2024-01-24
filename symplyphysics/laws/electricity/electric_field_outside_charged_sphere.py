from sympy import (Eq, solve, pi)
from sympy.physics.units import electric_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Let there be a uniformly charged sphere. Outside a uniformly charged sphere, the electric field
## is exactly the same as that created by a point charge placed in the center of the sphere, equal
## in magnitude to the total charge of the sphere. Inside a uniformly charged sphere, the electric
## field intensity is zero.

## Law is: E = q / (4 * pi * e0 * r^2), where
## E - electric field intensity,
## e0 - electric constant,
## q - electric charge,
## r - distance.

electric_intensity = Symbol("electric_intensity", units.voltage / units.length)

charge = Symbol("charge", units.charge)
distance = Symbol("distance", units.length)

law = Eq(electric_intensity, charge / (4 * pi * electric_constant * distance**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(charge_=charge, distance_=distance)
@validate_output(electric_intensity)
def calculate_electric_intensity(charge_: Quantity, distance_: Quantity) -> Quantity:
    result_expr = solve(law, electric_intensity, dict=True)[0][electric_intensity]
    result_expr = result_expr.subs({
        charge: charge_,
        distance: distance_,
    })
    return Quantity(result_expr)
