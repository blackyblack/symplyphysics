from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## The current density vector J is proportional to the negative of the gradient of the neutron flux.
## The proportionality constant is called the diffusion coefficient and is denoted by the symbol D.

# Conditions
## - Low absorbing medium (Σa << Σs)
##   Fick’s law derivation assumes that the neutron flux, φ, is slowly varying.
##   Large variations in φ occur when Σa (neutron absorption) is large (compared to Σs).

## Law: D = 1 / 3 * Σtr
## Where:
## Σtr (macroscopic transport cross-section) is the macroscopic transport cross-section.
##   See [macroscopic transport cross-section](./macroscopic_transport_cross_section.py) implementation.
## D is the diffusion coefficient.

macroscopic_transport_cross_section = symbols('macroscopic_transport_cross_section')
neutron_diffusion_coefficient = symbols('neutron_diffusion_coefficient')

law = Eq(neutron_diffusion_coefficient, 1 / (3 * macroscopic_transport_cross_section))

def print():
    return pretty(law, use_unicode=False)

@validate_input(macroscopic_transport_cross_section_=(1 / units.length))
@validate_output(units.length)
def calculate_diffusion_coefficient(macroscopic_transport_cross_section_: Quantity) -> Quantity:
    result_coefficient_expr = solve(law, neutron_diffusion_coefficient, dict=True)[0][neutron_diffusion_coefficient]
    result_expr = result_coefficient_expr.subs({
        macroscopic_transport_cross_section: macroscopic_transport_cross_section_})
    return expr_to_quantity(result_expr, 'neutron_diffusion_coefficient')
