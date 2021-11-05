from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

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

diffusion_coefficient = symbols('diffusion_coefficient')
macroscopic_absorption_cross_section = symbols('macroscopic_absorption_cross_section')
diffusion_area = symbols('diffusion_area')

law = Eq(diffusion_area, diffusion_coefficient / macroscopic_absorption_cross_section)

def print():
    return pretty(law, use_unicode=False)

@validate_input(diffusion_coefficient_=units.length, macroscopic_absorption_cross_section_=(1 / units.length))
@validate_output(units.length**2)
def calculate_diffusion_area(
    diffusion_coefficient_: Quantity,
    macroscopic_absorption_cross_section_: Quantity) -> Quantity:
    result_diffusion_expr = solve(law, diffusion_area, dict=True)[0][diffusion_area]
    result_expr = result_diffusion_expr.subs({
        diffusion_coefficient: diffusion_coefficient_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_})
    return expr_to_quantity(result_expr, 'neutron_diffusion_area')
