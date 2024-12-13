"""
Entropy is derivative of free energy
====================================

Entropy of a system can be found if its free energy is known as a function of temperature.

**Links:**

#. `Wikipedia, follows from the corresponding fundamental relation <https://en.wikipedia.org/wiki/Fundamental_thermodynamic_relation>`__.
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
from symplyphysics.laws.thermodynamics import free_energy_differential

entropy = symbols.entropy
"""
:symbols:`entropy` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

free_energy = clone_as_function(
    symbols.helmholtz_free_energy,
    [temperature, volume, particle_count],
)
"""
:symbols:`helmholtz_free_energy` of the system as a function of :attr:`~temperature`,
:attr:`~volume`, and :attr:`~particle_count`.
"""

law = Eq(entropy, -1 * Derivative(free_energy(temperature, volume, particle_count), temperature))
"""
:laws:symbol::

:laws:latex::
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
