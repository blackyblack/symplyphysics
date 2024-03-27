from sympy import Eq, solve, ln
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output,)

## Description
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.
## An electrical breakdown is a phenomenon of a sharp increase in current that occurs when the field intensity is higher
## than the critical one - dielectric breakdown intensity.


## Law is: Umax = E * D * ln(D / d) / 2, where
## Umax - maximum voltage between the central conductor and the outer conductor,
## E - dielectric breakdown intensity,
## D - diameter of the outer conductor,
## d - diameter of the inner conductor.


maximum_voltage = Symbol("maximum_voltage", units.voltage)

breakdown_intensity = Symbol("breakdown_intensity", units.voltage / units.length)
outer_diameter = Symbol("outer_diameter", units.length)
inner_diameter = Symbol("inner_diameter", units.length)

law = Eq(maximum_voltage, breakdown_intensity * outer_diameter * ln(outer_diameter / inner_diameter) / 2)


def print_law() -> str:
    return print_expression(law)


@validate_input(breakdown_intensity_=breakdown_intensity, outer_diameter_=outer_diameter, inner_diameter_=inner_diameter)
@validate_output(maximum_voltage)
def calculate_maximum_voltage(breakdown_intensity_: Quantity, outer_diameter_: Quantity, inner_diameter_: Quantity) -> Quantity:
    if outer_diameter_.scale_factor <= inner_diameter_.scale_factor:
        raise ValueError("The outer diameter must be greater than the inner diameter")
    result_velocity_expr = solve(law, maximum_voltage, dict=True)[0][maximum_voltage]
    result_expr = result_velocity_expr.subs({
        breakdown_intensity: breakdown_intensity_,
        outer_diameter: outer_diameter_,
        inner_diameter: inner_diameter_
    })
    return Quantity(result_expr)
