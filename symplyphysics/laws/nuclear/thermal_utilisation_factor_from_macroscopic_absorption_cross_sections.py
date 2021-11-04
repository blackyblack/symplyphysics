from sympy.core.singleton import S
from symplyphysics import (
    symbols, Eq, pretty, solve, Probability, Quantity, units,
    validate_input, expr_to_quantity, convert_to
)

# Description
## Thermal neutron utilization factor (f), is the ratio of the number of neutrons absorbed in the fuel
## versus the total number of absorptions everywhere in the reactor, i.e., in the fuel, moderator, cladding,
## and other reactor materials

## Law: f = Σa_fuel / Σa_total
## Where:
## Σa_fuel - macroscopic absorption cross-section of the fuel.
## Σa_total - macroscopic absorption cross-section of the fuel, moderator, cladding, etc, combined.
## f - thermal neutron utilisation factor

macroscopic_fuel_absorption_cross_section = symbols('macroscopic_fuel_absorption_cross_section')
macroscopic_total_absorption_cross_section = symbols('macroscopic_total_absorption_cross_section')
thermal_utilisation_factor = symbols('thermal_utilisation_factor')

law = Eq(thermal_utilisation_factor,
    macroscopic_fuel_absorption_cross_section / macroscopic_total_absorption_cross_section)

def print():
    return pretty(law, use_unicode=False)

@validate_input(macroscopic_fuel_absorption_cross_section_=(1 / units.length), macroscopic_total_absorption_cross_section_=(1 / units.length))
def calculate_utilisation_factor(
    macroscopic_fuel_absorption_cross_section_: Quantity,
    macroscopic_total_absorption_cross_section_: Quantity) -> Probability:
    if macroscopic_fuel_absorption_cross_section_.scale_factor > macroscopic_total_absorption_cross_section_.scale_factor:
        raise ValueError(f"macroscopic_fuel_absorption_cross_section_ ({macroscopic_fuel_absorption_cross_section_.scale_factor}) should be <= "
        f"macroscopic_total_absorption_cross_section_ ({macroscopic_total_absorption_cross_section_.scale_factor})")

    result_factor_expr = solve(law, thermal_utilisation_factor, dict=True)[0][thermal_utilisation_factor]
    result_expr = result_factor_expr.subs({
        macroscopic_fuel_absorption_cross_section: macroscopic_fuel_absorption_cross_section_,
        macroscopic_total_absorption_cross_section: macroscopic_total_absorption_cross_section_})
    result_factor = expr_to_quantity(result_expr, 'utilisation_factor')
    return Probability(convert_to(result_factor, S.One).n())
