from sympy.functions import exp
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units, S,
    Probability, validate_input, expr_to_quantity, convert_to
)

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

absorber_atomic_number_density = symbols('absorber_atomic_number_density')
effective_resonance_integral = symbols('effective_resonance_integral')
average_lethargy_change = symbols('average_lethargy_change')
macroscopic_scattering_cross_section_moderator = symbols('macroscopic_scattering_cross_section_moderator')
resonance_escape_probability = symbols('resonance_escape_probability')

law = Eq(resonance_escape_probability,
    exp(-1 * (absorber_atomic_number_density * effective_resonance_integral) /
    (average_lethargy_change * macroscopic_scattering_cross_section_moderator)))

def print():
    return pretty(law, use_unicode=False)

@validate_input(
    absorber_atomic_number_density_=(1 / units.length ** 3),
    effective_resonance_integral_=(units.length ** 2),
    macroscopic_scattering_cross_section_moderator_=(1 / units.length))
def calculate_resonance_escape_probability(
    absorber_atomic_number_density_: Quantity,
    effective_resonance_integral_: Quantity,
    average_lethargy_change_: float,
    macroscopic_scattering_cross_section_moderator_: Quantity) -> Probability:

    result_factor_expr = solve(law, resonance_escape_probability, dict=True)[0][resonance_escape_probability]
    result_expr = result_factor_expr.subs({
        absorber_atomic_number_density: absorber_atomic_number_density_,
        effective_resonance_integral: effective_resonance_integral_,
        average_lethargy_change: average_lethargy_change_,
        macroscopic_scattering_cross_section_moderator: macroscopic_scattering_cross_section_moderator_})
    result_factor = expr_to_quantity(result_expr, 'resonance_escape_factor')
    return Probability(convert_to(result_factor, S.One).n())
