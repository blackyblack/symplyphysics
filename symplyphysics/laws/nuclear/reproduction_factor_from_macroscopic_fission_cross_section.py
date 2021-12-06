from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units, S,
    validate_input, expr_to_quantity, convert_to
)

# Description
## The reproduction factor η represents the number of fast neutrons produced per thermal neutron absorbed in the fuel.

## Law: η = v * Σf / (Σf + Σa)
## Where:
## v - average number of neutrons produced per fission.
## Σf - macroscopic fission cross-section of the fuel.
##   See [macroscopic cross-section](./macroscopic_cross_section_from_free_mean_path.py) implementation.
## Σa - macroscopic absorption cross-section of the fuel.
## η - neutron reproduction factor

neutrons_per_fission = symbols('neutrons_per_fission')
macroscopic_fission_cross_section = symbols('macroscopic_fission_cross_section')
macroscopic_absorption_cross_section = symbols('macroscopic_absorption_cross_section')
neutron_reproduction_factor = symbols('neutron_reproduction_factor')

law = Eq(neutron_reproduction_factor,
    neutrons_per_fission * macroscopic_fission_cross_section / macroscopic_absorption_cross_section)

def print():
    return pretty(law, use_unicode=False)

@validate_input(macroscopic_fission_cross_section_=(1 / units.length), macroscopic_absorption_cross_section_=(1 / units.length))
def calculate_reproduction_factor(
    neutrons_per_fission_: float,
    macroscopic_fission_cross_section_: Quantity,
    macroscopic_absorption_cross_section_: Quantity) -> float:
    result_factor_expr = solve(law, neutron_reproduction_factor, dict=True)[0][neutron_reproduction_factor]
    result_expr = result_factor_expr.subs({
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fission_cross_section: macroscopic_fission_cross_section_,
        macroscopic_absorption_cross_section: macroscopic_absorption_cross_section_})
    result_factor = expr_to_quantity(result_expr, 'reproduction_factor')
    return convert_to(result_factor, S.One).n()
