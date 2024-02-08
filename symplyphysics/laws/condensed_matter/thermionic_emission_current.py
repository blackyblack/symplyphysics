from sympy import (Eq, solve, exp)
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
## Thermionic emission is the liberation of electrons from an electrode by virtue of its temperature.
## This occurs because the thermal energy given to the charge carrier overcomes the work function of the material.

## Law is: j = A * T^2 * exp(-F / kT), where
## j - the current density of thermionic emission,
## A - Richardson's constant,
## F - work function of the metal,
## k - Boltzmann constant,
## T - temperature.

density_current = Symbol("density_current", units.current / units.area)

thermodynamic_work = Symbol("thermodynamic_work", units.energy)
temperature = Symbol("temperature", units.temperature)

richardson_constant = Quantity(120.17 * (units.ampere / units.kelvin**2 / units.centimeter**2))

law = Eq(
    density_current,
    richardson_constant * (temperature**2) * exp(-thermodynamic_work / (boltzmann * temperature)))


def print_law() -> str:
    return print_expression(law)


@validate_input(thermodynamic_work_=thermodynamic_work, temperature_=temperature)
@validate_output(density_current)
def calculate_current(thermodynamic_work_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = solve(law, density_current, dict=True)[0][density_current]
    result_expr = result_expr.subs({
        thermodynamic_work: thermodynamic_work_,
        temperature: temperature_,
    })
    return Quantity(result_expr)
