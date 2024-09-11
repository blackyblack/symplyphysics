"""
Entropy is derivative of Gibbs energy
=====================================

Entropy of a system can be found if its Gibbs energy is known as a function of temperature.
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
from symplyphysics.laws.thermodynamics import gibbs_energy_differential

entropy = Symbol("entropy", units.energy / units.temperature)
"""
Entropy of the system.

Symbol:
    :code:`S`
"""

gibbs_energy = Function("gibbs_energy", units.energy)
"""
Gibbs energy of the system.

Symbol:
    :code:`G(T, p, N)`
"""

temperature = symbols.temperature
"""
:attr:`~symplyphysics.symbols.temperature` of the system.
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure inside the system.

Symbol:
    :code:`p`
"""

particle_count = Symbol("particle_count", dimensionless)
"""
Number of particles in the system.

Symbol:
    :code:`N`
"""

law = Eq(entropy, -1 * Derivative(gibbs_energy(temperature, pressure, particle_count), temperature))
r"""
:code:`S = -1 * Derivative(G(T, p, N), T)`

Latex:
    .. math::
        S = - \left( \frac{\partial G}{\partial T} \right)_{p, N}
"""

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
