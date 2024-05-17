from sympy import Eq, Derivative, solve
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import internal_energy_differential

# Description
## Temperature of a thermodynamic system can be found when internal energy is known as a function of entropy.

# Law: T = (dU/dS)_(V, N)
## T - absolute temperature
## U - internal energy
## S - entropy
## V - volume
## N - particle count
## (d/dS)_(V, N) - derivative with respect to entropy at constant volume and particle count

temperature = symbols.thermodynamics.temperature
internal_energy = Function("internal_energy", units.energy)
entropy = Symbol("entropy", units.energy / units.temperature)
volume = Symbol("volume", units.volume)
particle_count = Symbol("particle_count", dimensionless)

law = Eq(temperature, Derivative(internal_energy(entropy, volume, particle_count), entropy))

# Derive from fundamental relation for internal energy

_internal_energy_change_eqn = internal_energy_differential.law.subs({
    internal_energy_differential.volume_change: 0,
    internal_energy_differential.particle_count_change: 0,
})

_temperature_derived = solve(_internal_energy_change_eqn,
    internal_energy_differential.temperature)[0].subs(
    internal_energy_differential.internal_energy_change,
    Derivative(internal_energy(entropy, volume, particle_count), entropy) *
    internal_energy_differential.entropy_change,
    )

assert expr_equals(_temperature_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    internal_energy_change_=internal_energy,
    entropy_change_=entropy,
)
@validate_output(temperature)
def calculate_temperature(
    internal_energy_change_: Quantity,
    entropy_change_: Quantity,
) -> Quantity:
    internal_energy_ = (internal_energy_change_ / entropy_change_) * entropy
    result = law.rhs.subs(internal_energy(entropy, volume, particle_count), internal_energy_).doit()
    return Quantity(result)
