"""
Entropy is derivative of free energy
====================================

Entropy of a system can be found if its free energy is known as a function of temperature.
"""

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
from symplyphysics.laws.thermodynamics import free_energy_differential

entropy = Symbol("entropy", units.energy / units.temperature)
"""
Entropy of the system.

Symbol:
    :code:`S`
"""

free_energy = Function("free_energy", units.energy)
"""
Free energy of the system.

Symbol:
    :code:`F(T, V, N)`
"""

temperature = symbols.temperature
"""
:attr:`~symplyphysics.symbols.temperature` of the system.
"""

volume = Symbol("volume", units.volume)
"""
Volume of the system.

Symbol:
    :code:`V`
"""

particle_count = Symbol("particle_count", dimensionless)
"""
Number of particles in the system.

Symbol:
    :code:`N`
"""

law = Eq(entropy, -1 * Derivative(free_energy(temperature, volume, particle_count), temperature))
r"""
:code:`S = -1 * Derivative(F(T, V, N), T)`

Latex:
    .. math::
        S = - \left( \frac{\partial F}{\partial T} \right)_{V, N}
"""

# Derive from expression of free energy differential

_free_energy_change_eqn = free_energy_differential.law.subs({
    free_energy_differential.volume_change: 0,
    free_energy_differential.particle_count_change: 0,
})

_entropy_derived = solve(_free_energy_change_eqn, free_energy_differential.entropy)[0].subs(
    free_energy_differential.free_energy_change,
    Derivative(free_energy(temperature, volume, particle_count), temperature) *
    free_energy_differential.temperature_change)

assert expr_equals(_entropy_derived, law.rhs)


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
    result = law.rhs.subs(free_energy(temperature, volume, particle_count), free_energy_).doit()
    return Quantity(result)
