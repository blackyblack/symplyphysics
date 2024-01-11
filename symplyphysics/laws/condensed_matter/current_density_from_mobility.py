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
## Current density is the amount of charge per unit time that flows through a unit area of a chosen
## cross section. The current density vector is defined as a vector whose magnitude is the electric
## current per cross-sectional area at a given point in space, its direction being that of the motion
## of the positive charges at this point. In SI base units, the electric current density is measured
## in amperes per square metre.
## Mobility - characterises how quickly an electron or hole can move through a metal or semiconductor
## when pulled by an electric fiel.

## Law is: j = q * (-n1 * u1 + n2 * u2) * E, where
## j - current density of charge carriers,
## q - electron charge modulus,
## n1 - concentration of electrons,
## n2 - concentration of holes,
## u1 - mobility of electrons,
## u2 - mobility of holes,
## E - electric intensity (physical field that surrounds electrically charged particles).

density_current = Symbol("density_current", units.current / units.area)

electrons_concentration = Symbol("electrons_concentration", 1 / units.volume)
holes_concentration = Symbol("holes_concentration", 1 / units.volume)
electrons_mobility = Symbol("electrons_mobility", units.area / units.voltage / units.time)
holes_mobility = Symbol("holes_mobility", units.area / units.voltage / units.time)
electric_intensity = Symbol("electric_intensity", units.voltage / units.length)

charge = Quantity(1.6e-19 * units.coulomb)

law = Eq(
    density_current,
    charge *
    (-electrons_concentration * electrons_mobility + holes_concentration * holes_mobility) *
    electric_intensity)


def print_law() -> str:
    return print_expression(law)


@validate_input(electrons_concentration_=electrons_concentration,
    holes_concentration_=holes_concentration,
    electrons_mobility_=electrons_mobility,
    holes_mobility_=holes_mobility,
    electric_intensity_=electric_intensity)
@validate_output(density_current)
def calculate_current_density(electrons_concentration_: Quantity, holes_concentration_: Quantity,
    electrons_mobility_: Quantity, holes_mobility_: Quantity,
    electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, density_current, dict=True)[0][density_current]
    result_expr = result_expr.subs({
        electrons_concentration: electrons_concentration_,
        holes_concentration: holes_concentration_,
        electrons_mobility: electrons_mobility_,
        holes_mobility: holes_mobility_,
        electric_intensity: electric_intensity_
    })
    return Quantity(result_expr)
