from sympy import Eq, solve
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
## Thermal neutron utilization factor (f), is the ratio of the number of neutrons absorbed in the fuel
## versus the total number of absorptions everywhere in the reactor, i.e., in the fuel, moderator, cladding,
## and other reactor materials

## Law: f = Σa_fuel / Σa_total
## Where:
## Σa_fuel - macroscopic absorption cross-section of the fuel.
##   See [macroscopic cross-section](./macroscopic_cross_section_from_free_mean_path.py) implementation.
## Σa_total - macroscopic absorption cross-section of the fuel, moderator, cladding, etc, combined.
## f - thermal neutron utilisation factor

# Links:
## ScienceDirect <https://www.sciencedirect.com/topics/engineering/thermal-utilisation-factor>
## NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/nuclear-fission-chain-reaction/thermal-utilization-factor/>
## Wikipedia, second row in table <https://en.wikipedia.org/wiki/Six_factor_formula>

macroscopic_fuel_absorption_cross_section = Symbol("macroscopic_fuel_absorption_cross_section",
    1 / units.length)
macroscopic_total_absorption_cross_section = Symbol("macroscopic_total_absorption_cross_section",
    1 / units.length)
thermal_utilisation_factor = Symbol("thermal_utilisation_factor", dimensionless)

law = Eq(thermal_utilisation_factor,
    macroscopic_fuel_absorption_cross_section / macroscopic_total_absorption_cross_section)


@validate_input(
    macroscopic_fuel_absorption_cross_section_=macroscopic_fuel_absorption_cross_section,
    macroscopic_total_absorption_cross_section_=macroscopic_total_absorption_cross_section)
@validate_output(thermal_utilisation_factor)
def calculate_utilisation_factor(
        macroscopic_fuel_absorption_cross_section_: Quantity,
        macroscopic_total_absorption_cross_section_: Quantity) -> Probability:
    if macroscopic_fuel_absorption_cross_section_.scale_factor > macroscopic_total_absorption_cross_section_.scale_factor:
        raise ValueError(
            f"macroscopic_fuel_absorption_cross_section_ ({macroscopic_fuel_absorption_cross_section_.scale_factor}) should be <= "
            f"macroscopic_total_absorption_cross_section_ ({macroscopic_total_absorption_cross_section_.scale_factor})"
        )

    result_factor_expr = solve(law, thermal_utilisation_factor,
        dict=True)[0][thermal_utilisation_factor]
    result_expr = result_factor_expr.subs({
        macroscopic_fuel_absorption_cross_section: macroscopic_fuel_absorption_cross_section_,
        macroscopic_total_absorption_cross_section: macroscopic_total_absorption_cross_section_
    })
    return Probability(convert_to_float(result_expr))
