from sympy import (Eq, solve, exp)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless
)

# Description
## Thermionic emission is the liberation of electrons from an electrode by virtue of its temperature.
## This occurs because the thermal energy given to the charge carrier overcomes the work function of the material.

## Law is: j = A * T^2 * exp(-F / kT), where
## j - the current density of thermionic emission,
## A - Richardson's constant,
## F - work function of the metall,
## k - Boltzmann constant,
## T - temperature.

density_current = Symbol("density_current", units.current / units.length**2)

thermodynamic_work = Symbol("thermodynamic_work", units.energy)
electron_mass_ratio = Symbol("electron_mass_ratio", dimensionless)
temperature = Symbol("temperature", units.temperature)

richardson_constant = Quantity(120 * (units.A / units.kelvin**2 / units.centimeter**2))
boltzmann_constant = Quantity(1.380649e-23 * (units.joule / units.kelvin))

law = Eq(density_current, richardson_constant * electron_mass_ratio * (temperature**2) * exp(-thermodynamic_work / (boltzmann_constant * temperature)))


def print_law() -> str:
    return print_expression(law)


@validate_input(thermodynamic_work_=thermodynamic_work, electron_mass_ratio_=electron_mass_ratio, temperature_=temperature)
@validate_output(density_current)
def calculate_current(thermodynamic_work_: Quantity, electron_mass_ratio_: Quantity, temperature_: Quantity) -> Quantity:
    result_momentum_expr = solve(law, density_current, dict=True)[0][density_current]
    result_expr = result_momentum_expr.subs({
        thermodynamic_work: thermodynamic_work_,
        electron_mass_ratio: electron_mass_ratio_,
        temperature: temperature_,

    })
    return Quantity(result_expr)
