from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The current density vector J is proportional to the negative of the gradient of the neutron flux.
## The proportionality constant is called the diffusion coefficient and is denoted by the symbol D.

# Conditions
## - Low absorbing medium (Σa << Σs)
##   Fick’s law derivation assumes that the neutron flux, φ, is slowly varying.
##   Large variations in φ occur when Σa (macroscopic absorption cross-section) is large compared to
##   Σs (macroscopic scattering cross-section).

## Law: D = 1 / 3 * Σtr
## Where:
## Σtr (macroscopic transport cross-section) is the macroscopic transport cross-section.
##   See [macroscopic transport cross-section](./macroscopic_transport_cross_section.py) implementation.
## D is the diffusion coefficient.

# Links:
## NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/diffusion-coefficient/>

macroscopic_transport_cross_section = Symbol("macroscopic_transport_cross_section",
    1 / units.length)
neutron_diffusion_coefficient = Symbol("neutron_diffusion_coefficient", units.length)

law = Eq(neutron_diffusion_coefficient, 1 / (3 * macroscopic_transport_cross_section))


@validate_input(macroscopic_transport_cross_section_=macroscopic_transport_cross_section)
@validate_output(neutron_diffusion_coefficient)
def calculate_diffusion_coefficient(macroscopic_transport_cross_section_: Quantity) -> Quantity:
    result_coefficient_expr = solve(law, neutron_diffusion_coefficient,
        dict=True)[0][neutron_diffusion_coefficient]
    result_expr = result_coefficient_expr.subs(
        {macroscopic_transport_cross_section: macroscopic_transport_cross_section_})
    return Quantity(result_expr)
