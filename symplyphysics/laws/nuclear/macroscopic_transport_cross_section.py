from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    dimensionless,
    validate_input,
    validate_output,
)

# Description
## The transport mean free path (λtr) is an average distance a neutron will move in its original direction
## after infinite number of scattering collisions.

## Macroscopic transport cross-section: Σtr = Σs * (1 - μ)
## Where:
## Σs (macroscopic scattering cross-section) is the macroscopic cross-section for scattering interaction.
##   See [macroscopic cross-section](./macroscopic_cross_section_from_microscopic_cross_section.py) implementation.
## μ (average scattering angle cosine) is average value of the cosine of the angle in the lab system at which
##   neutrons are scattered in the medium.
##   See [average scattering angle cosine](./most_neutron_energies_scattering_angle_average_cosine.py) implementation.
## Σtr (macroscopic transport cross-section) is the macroscopic cross-section for transport mean free path.

# Links:
## NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/nuclear-engineering-fundamentals/neutron-nuclear-reactions/macroscopic-cross-section/>
## ScienceDirect <https://www.sciencedirect.com/topics/engineering/macroscopic-cross-section>
## NOTE find more fitting links

macroscopic_scattering_cross_section = Symbol("macroscopic_scattering_cross_section",
    1 / units.length)
average_scattering_angle_cosine = Symbol("average_scattering_angle_cosine", dimensionless)
macroscopic_transport_cross_section = Symbol("macroscopic_transport_cross_section",
    1 / units.length)

law = Eq(macroscopic_transport_cross_section,
    macroscopic_scattering_cross_section * (1 - average_scattering_angle_cosine))


@validate_input(macroscopic_scattering_cross_section_=macroscopic_scattering_cross_section,
    average_scattering_angle_cosine_=average_scattering_angle_cosine)
@validate_output(macroscopic_transport_cross_section)
def calculate_cross_section(macroscopic_scattering_cross_section_: Quantity,
    average_scattering_angle_cosine_: float) -> Quantity:
    result_cross_section_expr = solve(law, macroscopic_transport_cross_section,
        dict=True)[0][macroscopic_transport_cross_section]
    result_expr = result_cross_section_expr.subs({
        macroscopic_scattering_cross_section: macroscopic_scattering_cross_section_,
        average_scattering_angle_cosine: average_scattering_angle_cosine_
    })
    return Quantity(result_expr)
