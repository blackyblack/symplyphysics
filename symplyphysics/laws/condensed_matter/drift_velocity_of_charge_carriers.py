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
## Drift velocity is the average velocity attained by charged particles, such as electrons,
## in a material due to an electric field. In general, an electron in a conductor will propagate
## randomly at the Fermi velocity, resulting in an average velocity of zero. Applying an electric
## field adds to this random motion a small net flow in one direction; this is the drift.
## Mobility - characterises how quickly an electron or hole can move through a metal or semiconductor
## when pulled by an electric fiel.

## Law is: v = u * E, where
## v - drift velocity of charge carriers (electrons or holes),
## u - mobility of charge carriers (electrons or holes),
## E - electric field intensity.

drift_velocity = Symbol("drift_velocity", units.velocity)

charge_carriers_mobility = Symbol("charge_carriers_mobility", units.area / units.voltage / units.time)
electric_intensity = Symbol("electric_intensity", units.voltage / units.length)

law = Eq(drift_velocity, charge_carriers_mobility * electric_intensity)


def print_law() -> str:
    return print_expression(law)


@validate_input(charge_carriers_mobility_=charge_carriers_mobility, electric_intensity_=electric_intensity)
@validate_output(drift_velocity)
def calculate_velocity(charge_carriers_mobility_: Quantity, electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, drift_velocity, dict=True)[0][drift_velocity]
    result_expr = result_expr.subs({
        charge_carriers_mobility: charge_carriers_mobility_,
        electric_intensity: electric_intensity_,
    })
    return Quantity(result_expr)
