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
## The electric dipole moment is a measure of the separation of positive and negative electrical charges within a system,
## that is, a measure of the system's overall polarity. The SI unit for electric dipole moment is the coulomb-meter.
## The debye (D) is another unit of measurement used in atomic physics and chemistry.

## Law is: p = q * l, where
## p - electric dipole moment,
## q - electric charge,
## l - distance between charges.

electric_moment = Symbol("electric_moment", units.charge * units.length)

charge = Symbol("charge", units.charge)
distance = Symbol("distance", units.length)

law = Eq(
    electric_moment,
    charge * distance)


def print_law() -> str:
    return print_expression(law)


@validate_input(charge_=charge, distance_=distance)
@validate_output(electric_moment)
def calculate_electric_moment(charge_: Quantity, distance_: Quantity) -> Quantity:
    result_expr = solve(law, electric_moment, dict=True)[0][electric_moment]
    result_expr = result_expr.subs({
        charge: charge_,
        distance: distance_,
    })
    return Quantity(result_expr)
