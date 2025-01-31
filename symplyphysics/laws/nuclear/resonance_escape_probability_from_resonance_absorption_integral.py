from sympy import Eq, solve, exp
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    dimensionless,
    convert_to_float,
    validate_input,
    validate_output,
)
from symplyphysics.core.symbols.probability import Probability

# Description
## The resonance escape probability, symbolized by p, is the probability that a neutron will be
## slowed to thermal energy and escape resonance capture. This probability is defined as the ratio
## of the number of neutrons that reach thermal energies to the number of fast neutrons that slow down.

# Conditions
## - Homogeneous reactor.
## - Weak fast absorptions.
## - Very predominant absorber.

## Law: p ≈ e^(-(Na * Ieff) / (ξ * Σs_mod))
## Where:
## e - exponent.
## Na - atomic number density of the absorber.
##    See [atomic number density](./chemistry/atomic_number_density_from_material_density_atomic_weight.py) implementation.
## ξ - average lethargy gain per scattering event. Lethargy is a dimensionless logarithm of the ratio of the energy of
##    source neutrons to the energy of neutrons after a collision.
## Ieff - effective resonance integral. In practical situations, this integral strongly depends on the geometry of the unit cell.
##    It is too complex to calculate Ieff. Please refer to the tables and specific energy levels for proper data.
## Σs_mod - moderator's macroscopic scattering cross-section.
##   See [macroscopic cross-section](./macroscopic_cross_section_from_free_mean_path.py) implementation.
## p - resonance escape probability

# Links:
## Wikipedia, article <https://en.wikipedia.org/wiki/Resonance_escape_probability#Effective_resonance_integral>
## Wikipedia, third row in table <https://en.wikipedia.org/wiki/Six_factor_formula>

absorber_atomic_number_density = Symbol("absorber_atomic_number_density", 1 / units.volume)
effective_resonance_integral = Symbol("effective_resonance_integral", units.length**2)
average_lethargy_change = Symbol("average_lethargy_change", dimensionless)
macroscopic_scattering_cross_section_moderator = Symbol(
    "macroscopic_scattering_cross_section_moderator", 1 / units.length)
resonance_escape_probability = Symbol("resonance_escape_probability", dimensionless)

law = Eq(
    resonance_escape_probability,
    exp(-1 * (absorber_atomic_number_density * effective_resonance_integral) /
    (average_lethargy_change * macroscopic_scattering_cross_section_moderator)))


@validate_input(absorber_atomic_number_density_=absorber_atomic_number_density,
    effective_resonance_integral_=effective_resonance_integral,
    average_lethargy_change_=average_lethargy_change,
    macroscopic_scattering_cross_section_moderator_=macroscopic_scattering_cross_section_moderator)
@validate_output(resonance_escape_probability)
def calculate_resonance_escape_probability(
        absorber_atomic_number_density_: Quantity, effective_resonance_integral_: Quantity,
        average_lethargy_change_: float,
        macroscopic_scattering_cross_section_moderator_: Quantity) -> Probability:

    result_factor_expr = solve(law, resonance_escape_probability,
        dict=True)[0][resonance_escape_probability]
    result_expr = result_factor_expr.subs({
        absorber_atomic_number_density:
            absorber_atomic_number_density_,
        effective_resonance_integral:
            effective_resonance_integral_,
        average_lethargy_change:
            average_lethargy_change_,
        macroscopic_scattering_cross_section_moderator:
            macroscopic_scattering_cross_section_moderator_
    })
    return Probability(convert_to_float(result_expr))
