from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

# Description
## Quality factor is the property of oscillatiing system. It shows the ratio between amount of energy stored in system and power losses.

# Definition: Q = w * W / P.

# Where:
## Q is quality factor,
## w is raisonant circular frequency,
## W is energy stored in system,
## P is dissipated power.

quality_factor = Symbol("quality_factor", dimensionless)
raisonant_frequency = Symbol("raisonant_frequency", units.frequency)
stored_energy = Symbol("stored_energy", units.energy)
dissipated_power = Symbol("dissipated_power", units.power)

definition = Eq(quality_factor, raisonant_frequency * stored_energy / dissipated_power)

definition_units_SI = dimensionless

def print_law() -> str:
    return print_expression(definition)


@validate_input(frequency_=raisonant_frequency, energy_=stored_energy, power_=dissipated_power)
@validate_output(quality_factor)
def calculate_quality_factor(frequency_: Quantity, energy_: Quantity, power_: Quantity) -> Quantity:
    result_factor_expr = solve(definition, quality_factor, dict=True)[0][quality_factor]
    result_expr = result_factor_expr.subs({raisonant_frequency: frequency_, stored_energy: energy_, dissipated_power: power_})
    return Quantity(result_expr)
