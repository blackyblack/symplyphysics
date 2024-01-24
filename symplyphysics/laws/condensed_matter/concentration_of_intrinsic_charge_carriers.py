from sympy import (Eq, solve, sqrt, exp)
from sympy.physics.units import boltzmann
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## In the absence of external influences (lighting, electric field, etc.), there is a
## concentration of free charge carriers in the semiconductor.

## Law is: n = sqrt(Nc*Nv) * exp(-Eg / 2kT), where
## n is concentration of intrinsic charge carriers,
## Nc - effective density of states in the conduction band,
## Nv - effective density of states in the valence band,
## Eg - the width of band gap,
## k - Boltzmann constant,
## T - temperature.

charge_carriers_concentration = Symbol("charge_carriers_concentration", 1 / units.volume)

density_of_states_in_conduction_band = Symbol("density_of_states_in_conduction_band",
    1 / units.volume)
density_of_states_in_valence_band = Symbol("density_of_states_in_valence_band", 1 / units.volume)
temperature = Symbol("temperature", units.temperature)
band_gap = Symbol("band_gap", units.energy)

law = Eq(
    charge_carriers_concentration,
    sqrt(density_of_states_in_conduction_band * density_of_states_in_valence_band) * exp(-band_gap /
    (2 * boltzmann * temperature)))


def print_law() -> str:
    return print_expression(law)


@validate_input(density_of_states_in_conduction_band_=density_of_states_in_conduction_band,
    density_of_states_in_valence_band_=density_of_states_in_valence_band,
    band_gap_=band_gap,
    temperature_=temperature)
@validate_output(charge_carriers_concentration)
def calculate_concentration(density_of_states_in_conduction_band_: Quantity,
    density_of_states_in_valence_band_: Quantity, band_gap_: Quantity,
    temperature_: Quantity) -> Quantity:
    result_expr = solve(law, charge_carriers_concentration,
        dict=True)[0][charge_carriers_concentration]
    result_expr = result_expr.subs({
        density_of_states_in_conduction_band: density_of_states_in_conduction_band_,
        density_of_states_in_valence_band: density_of_states_in_valence_band_,
        band_gap: band_gap_,
        temperature: temperature_,
    })
    return Quantity(result_expr)
