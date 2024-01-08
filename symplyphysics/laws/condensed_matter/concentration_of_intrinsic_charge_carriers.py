from sympy import (Eq, solve, pi, sqrt, exp)
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

## Law is: n = ((Nc*Nv)^1/2) * exp(-Eg / 2kT), where
## n is concentration intrinsic charge carriers,
## Nc and Nv - effective densities of states in the conduction band and in the valence band,
## Eg - the width of band gap,
## k - Boltzmann constant,
## T - temperature.

concentration_intrinsic = Symbol("concentration_intrinsic", 1 / units.length**3)

density_conductivity = Symbol("density_conductivity", 1 / units.length**3)
density_valence = Symbol("density_valence", 1 / units.length**3)
temperature = Symbol("temperature", units.temperature)
band_gap = Symbol("band_gap", units.energy)

boltzmann_constant = Quantity(1.380649e-23 * (units.joule / units.kelvin))

law = Eq(concentration_intrinsic, sqrt(density_conductivity * density_valence) * exp(-band_gap/(2 * boltzmann_constant * temperature)))


def print_law() -> str:
    return print_expression(law)


@validate_input(density_conductivity_=density_conductivity, density_valence_=density_valence, band_gap_=band_gap, temperature_=temperature)
@validate_output(concentration_intrinsic)
def calculate_concentration(density_conductivity_: Quantity, density_valence_: Quantity, band_gap_: Quantity, temperature_: Quantity) -> Quantity:
    result_momentum_expr = solve(law, concentration_intrinsic, dict=True)[0][concentration_intrinsic]
    result_expr = result_momentum_expr.subs({
        density_conductivity: density_conductivity_,
        density_valence: density_valence_,
        band_gap: band_gap_,
        temperature: temperature_,

    })
    return Quantity(result_expr)
