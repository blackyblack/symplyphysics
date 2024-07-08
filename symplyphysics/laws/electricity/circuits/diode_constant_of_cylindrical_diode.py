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
## In a cylindrical diode, the cathode is located in the center, and the anode is located around it in the form of a cylinder.
## The radius of the anode is usually much larger than the radius of the cathode. The diode constant depends on the radii and on the area of the anode.

## Law is: g = (4 / 9) * e0 * sqrt(2 * e / m) * sa / (ra^2 * (1 - rc / ra)^2), where
## g - diode constant,
## e0 - electric constant,
## e - elementary charge,
## m - electron rest mass,
## sa - anode area,
## ra - anode radius,
## rc - cathode radius.

diode_constant = Symbol("diode_constant", units.current / units.voltage**Rational(3, 2))

anode_area = Symbol("anode_area", units.area)
anode_radius = Symbol("anode_radius", units.length)
cathode_radius = Symbol("cathode_radius", units.length)

law = Eq(diode_constant, Rational(4, 9) * electric_constant * sqrt(2 * elementary_charge / electron_rest_mass) * anode_area / (anode_radius**2 * (1 - cathode_radius / anode_radius)**2))


@validate_input(anode_area_=anode_area,
    anode_radius_=anode_radius,
    cathode_radius_=cathode_radius)
@validate_output(diode_constant)
def calculate_diode_constant(anode_area_: Quantity, anode_radius_: Quantity, cathode_radius_: Quantity) -> Quantity:
    if anode_radius_.scale_factor <= cathode_radius_.scale_factor:
        raise ValueError("The anode radius must be greater than the cathode radius")
    result_expr = solve(law, diode_constant,
        dict=True)[0][diode_constant]
    result_expr = result_expr.subs({
        anode_area: anode_area_,
        anode_radius: anode_radius_,
        cathode_radius: cathode_radius_,
    })
    return Quantity(result_expr)
