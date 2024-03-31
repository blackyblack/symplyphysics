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
## Intensity is a scalar physical quantity that quantitatively characterizes the power carried by
## a wave in the direction of propagation

## Law: I = E / (S * t), where
## I - intensity of incident radiation,
## E - energy of incident radiation,
## S - area,
## t - time.

intensity = Symbol("intensity", units.power / units.area)

energy = Symbol("energy", units.energy)
area = Symbol("area", units.area)
time = Symbol("time", units.time)

law = Eq(intensity, energy / (area * time))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    energy_=energy,
    area_=area,
    time_=time,
)
@validate_output(intensity)
def calculate_intensity(energy_: Quantity, area_: Quantity, time_: Quantity) -> Quantity:
    result_expr = solve(law, intensity, dict=True)[0][intensity]
    intensity_applied = result_expr.subs({energy: energy_, area: area_, time: time_})
    return Quantity(intensity_applied)
