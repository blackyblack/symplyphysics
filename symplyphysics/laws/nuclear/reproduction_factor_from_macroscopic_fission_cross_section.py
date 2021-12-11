from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units, S,
    validate_input, expr_to_quantity, convert_to
)

# Description
## The reproduction factor η represents the number of fast neutrons produced per thermal neutron absorbed in the fuel.

## Law: η = v * Σf_fuel / Σa_fuel
## Where:
## v - average number of neutrons produced per fission.
## Σf_fuel - macroscopic fission cross-section of the fuel.
##   See [macroscopic cross-section](./macroscopic_cross_section_from_free_mean_path.py) implementation.
## Σa_fuel - macroscopic absorption cross-section of the fuel. Equals to macroscopic fission (Σf_fuel) + macroscopic capture
##   (Σc_fuel) cross-sections.
## η - neutron reproduction factor

neutrons_per_fission = symbols('neutrons_per_fission')
macroscopic_fuel_fission_cross_section = symbols('macroscopic_fuel_fission_cross_section')
macroscopic_fuel_absorption_cross_section = symbols('macroscopic_fuel_absorption_cross_section')
neutron_reproduction_factor = symbols('neutron_reproduction_factor')

law = Eq(neutron_reproduction_factor,
    neutrons_per_fission * macroscopic_fuel_fission_cross_section / macroscopic_fuel_absorption_cross_section)

def print():
    return pretty(law, use_unicode=False)

@validate_input(macroscopic_fuel_fission_cross_section_=(1 / units.length), macroscopic_fuel_absorption_cross_section_=(1 / units.length))
def calculate_reproduction_factor(
    neutrons_per_fission_: float,
    macroscopic_fuel_fission_cross_section_: Quantity,
    macroscopic_fuel_absorption_cross_section_: Quantity) -> float:

    if macroscopic_fuel_fission_cross_section_.scale_factor > macroscopic_fuel_absorption_cross_section_.scale_factor:
        raise ValueError(f"macroscopic_fuel_fission_cross_section_ ({macroscopic_fuel_fission_cross_section_.scale_factor}) should be <= "
        f"macroscopic_fuel_absorption_cross_section_ ({macroscopic_fuel_absorption_cross_section_.scale_factor})")

    result_factor_expr = solve(law, neutron_reproduction_factor, dict=True)[0][neutron_reproduction_factor]
    result_expr = result_factor_expr.subs({
        neutrons_per_fission: neutrons_per_fission_,
        macroscopic_fuel_fission_cross_section: macroscopic_fuel_fission_cross_section_,
        macroscopic_fuel_absorption_cross_section: macroscopic_fuel_absorption_cross_section_})
    result_factor = expr_to_quantity(result_expr, 'reproduction_factor')
    return convert_to(result_factor, S.One).n()
