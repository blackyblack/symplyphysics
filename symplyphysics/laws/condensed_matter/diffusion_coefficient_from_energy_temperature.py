from sympy import (Eq, solve, exp)
from sympy.physics.units import boltzmann_constant
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The diffusion coefficient is a quantitative characteristic of the diffusion rate,
## equal to the amount of matter passing per unit time through a section of a unit area as
## a result of the thermal motion of molecules.

## Law is: D = D0 * exp(-Ea / (k * T)), where
## R - diffusion coefficient,
## D0 - diffusion constant (tabular value depends on dopant),
## Ea - activation energy of dopant,
## k - boltzmann constant,
## T - temperature.

# Conditions:
## 1. The material is a solid.

# Links: Wikipedia <https://en.wikipedia.org/wiki/Mass_diffusivity#Solids>

diffusion_coefficient = Symbol("diffusion_coefficient", units.area / units.time)

energy = Symbol("energy", units.energy)
diffusion_constant = Symbol("diffusion_constant", units.area / units.time)
temperature = symbols.temperature

law = Eq(diffusion_coefficient,
    diffusion_constant * exp(-energy / (boltzmann_constant * temperature)))


@validate_input(energy_=energy, diffusion_constant_=diffusion_constant, temperature_=temperature)
@validate_output(diffusion_coefficient)
def calculate_diffusion_coefficient(energy_: Quantity, diffusion_constant_: Quantity,
    temperature_: Quantity) -> Quantity:
    result_expr = solve(law, diffusion_coefficient, dict=True)[0][diffusion_coefficient]
    result_expr = result_expr.subs({
        energy: energy_,
        diffusion_constant: diffusion_constant_,
        temperature: temperature_
    })
    return Quantity(result_expr)
