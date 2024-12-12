"""
Temperature is derivative of internal energy
============================================

Temperature of a thermodynamic system can be found when internal energy is known as a function of entropy.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Internal_energy#Internal_energy_of_multi-component_systems>`__.
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
from symplyphysics.laws.thermodynamics import internal_energy_differential

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

entropy = symbols.entropy
"""
:symbols:`entropy` of the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

internal_energy = clone_as_function(
    symbols.internal_energy,
    [entropy, volume, particle_count],
)
"""
:symbols:`internal_energy` of the system as a function of its natural variables.
"""

law = Eq(temperature, Derivative(internal_energy(entropy, volume, particle_count), entropy))
"""
:laws:symbol::

:laws:latex::
"""

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
