from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)

# Description
## The physical meaning of the diffusion length can be seen by calculating the mean square distance that
## a neutron travels in the one direction from the plane source to its absorption point.

## Law: L^2 = D / Σa
## Where:
## Σa - macroscopic absorption cross-section of the fuel.
##   See [macroscopic transport cross-section](./macroscopic_transport_cross_section.py) implementation.
## D - diffusion coefficient.
##   See [diffusion coefficient](./neutron_diffusion_coefficient_from_scattering_cross_section.py) implementation.
## L^2 - diffusion area.
## L - diffusion length.

# Links:
## NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/diffusion-length/>

diffusion_coefficient = Symbol("diffusion_coefficient", units.length)
macroscopic_absorption_cross_section = Symbol("macroscopic_absorption_cross_section",
    1 / units.length)
diffusion_area = Symbol("diffusion_area", units.area)

law = Eq(diffusion_area, diffusion_coefficient / macroscopic_absorption_cross_section)


@validate_input(diffusion_coefficient_=diffusion_coefficient,
    macroscopic_absorption_cross_section_=macroscopic_absorption_cross_section)
@validate_output(diffusion_area)
def calculate_diffusion_area(diffusion_coefficient_: Quantity,
    macroscopic_absorption_cross_section_: Quantity) -> Quantity:
    result_diffusion_expr = solve(law, diffusion_area, dict=True)[0][diffusion_area]
    result_expr = result_diffusion_expr.subs({
        diffusion_coefficient: diffusion_coefficient_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_
    })
    return Quantity(result_expr)
