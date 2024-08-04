from sympy import Eq, Rational, solve, sqrt
from sympy.physics.units import electric_constant, elementary_charge, electron_rest_mass
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The current-voltage characteristic of a vacuum diode is described by the 3/2-power law. The diode
## constant in this law depends only on the relative position, shape and size of the electrodes of the vacuum diode.
## The vacuum diode may have different electrode geometries. Here we are talking about a vacuum diode with plane-parallel electrodes.

## Law is: g = (4 / 9) * e0 * sqrt(2 * e / m) * (s / d^2), where
## g - diode constant,
## e0 - electric constant,
## e - elementary charge,
## m - electron rest mass,
## s - area of one electrode,
## d - distance between electrodes.

diode_constant = Symbol("diode_constant", units.current / units.voltage**Rational(3, 2))

electrode_area = Symbol("electrode_area", units.area)
distance_between_electrodes = Symbol("distance_between_electrodes", units.length)

law = Eq(
    diode_constant,
    Rational(4, 9) * electric_constant * sqrt(2 * elementary_charge / electron_rest_mass) *
    (electrode_area / distance_between_electrodes**2))


@validate_input(electrode_area_=electrode_area,
    distance_between_electrodes_=distance_between_electrodes)
@validate_output(diode_constant)
def calculate_diode_constant(electrode_area_: Quantity,
    distance_between_electrodes_: Quantity) -> Quantity:
    result_expr = solve(law, diode_constant, dict=True)[0][diode_constant]
    result_expr = result_expr.subs({
        electrode_area: electrode_area_,
        distance_between_electrodes: distance_between_electrodes_,
    })
    return Quantity(result_expr)
