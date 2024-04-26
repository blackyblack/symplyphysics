from sympy import Eq, exp
from sympy.physics.units import boltzmann_constant, molar_gas_constant, planck
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

# Description
## To derive the heat capacity of a solid, one should account for quantum effects. Albert Einstein used the
## same model as in the classical case, namely the atoms are harmonic oscillators with three degrees
## of freedom, located in the nodes of the crystal lattice, performing thermal oscillations around the
## equlibrium positions with the same frequency. But he used a more correct expression for the energy
## of the oscillators, and although the result still only qualitatively describes the heat capacity of
## solids, it is a big achievement and the result has correct asymptotic behaviour for `T -> 0`.

## Law: C_V = 3 * R * x**2 * exp(x) / (exp(x) - 1)**2
## C_V - isochoric molar heat capacity
## R - molar gas constant
## h - Planck constant
## nu - frequency of thermal oscillations
## k - Boltzmann constant
## T - absolute temperature
## x = h * nu / (k * T) - ratio between photon and thermal energy

isochoric_molar_heat_capacity = Symbol(
    "isochoric_molar_heat_capacity",
    units.energy / (units.temperature * units.amount_of_substance)
)
frequency = Symbol("frequency", units.frequency)
temperature = symbols.thermodynamics.temperature

energy_ratio = (planck * frequency) / (boltzmann_constant * temperature)

law = Eq(
    isochoric_molar_heat_capacity,
    (3 * molar_gas_constant) * energy_ratio**2 * exp(energy_ratio) / (exp(energy_ratio) - 1)**2
)


@validate_input(
    frequency_=frequency,
    temperature_=temperature,
)
@validate_output(isochoric_molar_heat_capacity)
def calculate_isochoric_molar_heat_capacity(
    frequency_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        frequency: frequency_,
        temperature: temperature_,
    })
    return Quantity(result)
