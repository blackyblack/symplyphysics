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
## Entropy of a system can be found if its Gibbs energy is known as a function of temperature.

# Law: S = -(dG/dT)_(p, N)
## S - entropy
## G - Gibbs energy
## T - absolute temperature
## p - pressure
## N - particle count
## (d/dT)_(p, N) - derivative with respect to temperature at constant pressure and particle count

entropy = Symbol("entropy", units.energy / units.temperature)
gibbs_energy = Function("gibbs_energy", units.energy)
temperature = symbols.thermodynamics.temperature
pressure = Symbol("pressure", units.pressure)
particle_count = Symbol("particle_count", dimensionless)

law = Eq(entropy, -1 * Derivative(gibbs_energy(temperature, pressure, particle_count), temperature))


@validate_input(
    gibbs_energy_change_=gibbs_energy,
    temperature_change_=temperature,
)
@validate_output(entropy)
def calculate_entropy(
    gibbs_energy_change_: Quantity,
    temperature_change_: Quantity,
) -> Quantity:
    gibbs_energy_ = (gibbs_energy_change_ / temperature_change_) * temperature
    result = law.rhs.subs(
        gibbs_energy(temperature, pressure, particle_count), gibbs_energy_
    ).doit()
    return Quantity(result)
