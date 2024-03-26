from sympy import (Eq, solve, sqrt, ln,)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless,)


# Description
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.

## Law is: P = (U^2 / 120) *  sqrt(er / (mur * ln(D / d))), where
## P - power carried by the waveguide,
## U - voltage,
## er - relative permittivity of insulating material,
## mur - relative permeability of the insulator material,
## D - diameter of the outer conductor,
## d - diameter of the inner conductor.

waveguide_power = Symbol("waveguide_power", units.power)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
relative_permeability = Symbol("relative_permeability", dimensionless)
voltage = Symbol("voltage", units.voltage)
outer_diameter = Symbol("outer_diameter", units.length)
inner_diameter = Symbol("inner_diameter", units.length)

resistance = Quantity(120 * units.ohm)

law = Eq(waveguide_power, (voltage**2 / resistance) * sqrt(relative_permittivity / (relative_permeability * ln(outer_diameter / inner_diameter))))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability,
    voltage_=voltage,
    outer_diameter_=outer_diameter,
    inner_diameter_=inner_diameter)
@validate_output(waveguide_power)
def calculate_waveguide_power(relative_permittivity_: float, relative_permeability_: Quantity,
    voltage_: Quantity, outer_diameter_: Quantity,
    inner_diameter_: Quantity) -> Quantity:
    if outer_diameter_.scale_factor <= inner_diameter_.scale_factor:
        raise ValueError("The outer diameter must be greater than the inner diameter")
    result_expr = solve(law, waveguide_power, dict=True)[0][waveguide_power]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
        voltage: voltage_,
        outer_diameter: outer_diameter_,
        inner_diameter: inner_diameter_
    })
    return Quantity(result_expr)
