from sympy import Eq, Derivative, solve
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
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import gibbs_energy_differential

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

# Derive from the expression of Gibbs energy differential

_gibbs_energy_change_eqn = gibbs_energy_differential.law.subs({
    gibbs_energy_differential.pressure_change: 0,
    gibbs_energy_differential.particle_count_change: 0,
})

_entropy_derived = solve(
    _gibbs_energy_change_eqn,
    gibbs_energy_differential.entropy,
)[0].subs(
    gibbs_energy_differential.gibbs_energy_change,
    Derivative(gibbs_energy(temperature, pressure, particle_count), temperature) *
    gibbs_energy_differential.temperature_change)

assert expr_equals(_entropy_derived, law.rhs)


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
    result = law.rhs.subs(gibbs_energy(temperature, pressure, particle_count), gibbs_energy_).doit()
    return Quantity(result)
