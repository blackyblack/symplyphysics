from sympy import (Eq, solve, pi)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Hagen–Poiseuille law  is a physical law that gives the pressure drop in an incompressible and Newtonian
## fluid in laminar flow flowing through a long cylindrical pipe of constant cross section.

## Law: Q = (ΔP * pi * r^4) / 8 * eta * l
## Where:
## Q is flow of liquid
## ΔP is pressure differential between the two ends of the tube
## r is radius of the narrow tube
## eta is viscosity of liquid
## l is length of the arrow tube


flow_of_liquid = Symbol("flow_of_liquid", units.volume / units.time)
pressure_differential = Symbol("pressure_differential", units.pressure)
radius = Symbol("radius", units.length)
viscosity = Symbol("viscosity", units.pressure * units.time)
length = Symbol("length", units.length)

law = Eq(flow_of_liquid, (pressure_differential * pi * radius**4) / (8 * viscosity * length))


def print_law() -> str:
    return print_expression(law)


@validate_input(pressure_differential_=pressure_differential, radius_=radius, viscosity_=viscosity, length_=length)
@validate_output(flow_of_liquid)
def calculate_flow_of_liquid(pressure_differential_: Quantity, radius_: Quantity,
    viscosity_: Quantity, length_: Quantity) -> Quantity:
    result_expr = solve(law, flow_of_liquid, dict=True)[0][flow_of_liquid]
    result_flow_of_liquid = result_expr.subs({
        pressure_differential: pressure_differential_,
        radius: radius_,
        viscosity: viscosity_,
        length: length_,
    })
    return Quantity(result_flow_of_liquid)
