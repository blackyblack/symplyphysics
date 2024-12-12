from sympy import Eq, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## Consider a two-electrode gas-filled gap with flat electrodes. Let one of the electrodes emit charged particles into a gaseous medium.
## If the potential collecting these charged particles is set to another electrode, then a current will flow through the gap.

## Law is: E = (3 / 2) * sqrt(x / L) * U / L, where
## E - the electric intensity at the point between the electrodes,
## x - the coordinate of the point between the electrodes, located on the axis from the cathode to the anode,
## L - the distance between the electrodes,
## U - voltage between the electrodes.

# Conditions:
## - Emittivity of the emitter is unlimited (emission of charged particles is much greater than current between electrodes)
## - Gaseous medium between electrodes

# TODO: find link

electric_intensity = Symbol("electric_intensity", units.voltage / units.length)

coordinate = Symbol("coordinate", units.length)
distance_between_electrodes = Symbol("distance_between_electrodes", units.length)
voltage_between_electrodes = Symbol("voltage_between_electrodes", units.voltage)

law = Eq(electric_intensity, (3 / 2) * sqrt(coordinate / distance_between_electrodes) *
    voltage_between_electrodes / distance_between_electrodes)


@validate_input(coordinate_=coordinate,
    distance_between_electrodes_=distance_between_electrodes,
    voltage_between_electrodes_=voltage_between_electrodes)
@validate_output(electric_intensity)
def calculate_electric_intensity(coordinate_: Quantity, distance_between_electrodes_: Quantity,
    voltage_between_electrodes_: Quantity) -> Quantity:
    if coordinate_.scale_factor > distance_between_electrodes_.scale_factor:
        raise ValueError("The point of interest should be between electrodes")
    result_expr = solve(law, electric_intensity, dict=True)[0][electric_intensity]
    result_expr = result_expr.subs({
        coordinate: coordinate_,
        distance_between_electrodes: distance_between_electrodes_,
        voltage_between_electrodes: voltage_between_electrodes_,
    })
    return Quantity(result_expr)
