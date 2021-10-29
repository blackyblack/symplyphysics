from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## The transport mean free path (λtr) is an average distance a neutron will move in its original direction
## after infinite number of scattering collisions.

## Macroscopic transport cross-section: Σtr = Σs * (1 - μ)
## Where:
## Σs (macroscopic scattering cross-section) is the macroscopic cross-section for scattering interaction.
## μ (average_scattering_angle_cosine) is average value of the cosine of the angle in the lab system at which
##   neutrons are scattered in the medium.
##   See [most_neutron_energies_scattering_angle_average_cosine](./most_neutron_energies_scattering_angle_average_cosine.py)
## Σtr (macroscopic transport cross-section) is the macroscopic cross-section for transport mean free path.

macroscopic_scattering_cross_section = symbols('macroscopic_scattering_cross_section')
average_scattering_angle_cosine = symbols('average_scattering_angle_cosine')
macroscopic_transport_cross_section = symbols('macroscopic_transport_cross_section')

law = Eq(macroscopic_transport_cross_section,
    macroscopic_scattering_cross_section * (1 - average_scattering_angle_cosine))

def print():
    return pretty(law, use_unicode=False)

@validate_input(macroscopic_scattering_cross_section_=(1 / units.length))
@validate_output(1 / units.length)
def calculate_cross_section(macroscopic_scattering_cross_section_: Quantity, average_scattering_angle_cosine_: float) -> Quantity:
    result_cross_section_expr = solve(law, macroscopic_transport_cross_section, dict=True)[0][macroscopic_transport_cross_section]
    result_expr = result_cross_section_expr.subs({
        macroscopic_scattering_cross_section: macroscopic_scattering_cross_section_,
        average_scattering_angle_cosine: average_scattering_angle_cosine_})
    return expr_to_quantity(result_expr, 'macro_transport_cross_section')
