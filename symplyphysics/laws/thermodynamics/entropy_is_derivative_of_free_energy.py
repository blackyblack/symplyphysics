from sympy import Eq, Derivative
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
    symbols,
)

# Description
## Entropy of a system can be found if its free energy is known as a function of temperature.

# Law: S = -(dF/dT)_(V, N)
## S - entropy
## F - Helmholtz free energy
## T - absolute temperature
## V - volume
## N - particle count
## (d/dT)_(V, N) - derivative with respect to temperature at constant volume and particle count

entropy = Symbol("entropy", units.energy / units.temperature)
free_energy = Function("free_energy", units.energy)
temperature = symbols.thermodynamics.temperature
volume = Symbol("volume", units.volume)
particle_count = Symbol("particle_count", dimensionless)

law = Eq(entropy, -1 * Derivative(free_energy(temperature, volume, particle_count), temperature))


@validate_input(
    free_energy_change_=free_energy,
    temperature_change_=temperature,
)
@validate_output(entropy)
def calculate_entropy(
    free_energy_change_: Quantity,
    temperature_change_: Quantity,
) -> Quantity:
    free_energy_ = (free_energy_change_ / temperature_change_) * temperature
    result = law.rhs.subs(
        free_energy(temperature, volume, particle_count), free_energy_
    ).doit()
    return Quantity(result)
