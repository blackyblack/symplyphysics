"""
Entropy is derivative of Gibbs energy
=====================================

Entropy of a system can be found if its Gibbs energy is known as a function of temperature.
"""

from sympy import Eq, Derivative, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import gibbs_energy_differential

entropy = symbols.entropy
"""
:symbols:`entropy` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` inside the system.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

gibbs_energy = clone_as_function(
    symbols.gibbs_energy,
    [temperature, pressure, particle_count],
)
"""
:symbols:`gibbs_energy` of the system as a function of :attr:`~temperature`, :attr:`~pressure`,
and :attr:`~particle_count`.
"""

law = Eq(entropy, -1 * Derivative(gibbs_energy(temperature, pressure, particle_count), temperature))
"""
:laws:symbol::

:laws:latex::
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
