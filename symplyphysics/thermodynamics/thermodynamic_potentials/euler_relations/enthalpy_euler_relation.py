"""
Euler relation for enthalpy
===========================

The enthalpy differential can be integrated using the Euler's theorem on homogeneous functions for
internal energy.

**Notes:**

#. This formula works for a single-component system. For a multi-component system replace the
   product of chemical potential and particle count with a sum over each type of components.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Thermodynamic_potential#Euler_relations>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.thermodynamics.thermodynamic_potentials.definitions import enthalpy_is_internal_energy_plus_pressure_energy as enthalpy_def
from symplyphysics.thermodynamics.thermodynamic_potentials.euler_relations import internal_energy_euler_relation

enthalpy = symbols.enthalpy
"""
:symbols:`enthalpy` of the system. Also see :ref:`Enthalpy is internal energy plus pressure energy`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

entropy = symbols.entropy
"""
:symbols:`entropy` of the system. Also see :ref:`Entropy change in reversible process`.
"""

chemical_potential = symbols.chemical_potential
"""
:symbols:`chemical_potential` of the system. Also see
:ref:`Chemical potential is particle count derivative of internal energy`.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

law = Eq(enthalpy, temperature * entropy + chemical_potential * particle_count)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from enthalpy definition and internal energy Euler relation

_internal_energy_expr = internal_energy_euler_relation.law.rhs.subs({
    internal_energy_euler_relation.temperature: temperature,
    internal_energy_euler_relation.entropy: entropy,
    internal_energy_euler_relation.pressure: enthalpy_def.pressure,
    internal_energy_euler_relation.volume: enthalpy_def.volume,
    internal_energy_euler_relation.chemical_potential: chemical_potential,
    internal_energy_euler_relation.particle_count: particle_count,
})

_enthalpy_expr = enthalpy_def.law.rhs.subs({enthalpy_def.internal_energy: _internal_energy_expr})

assert expr_equals(_enthalpy_expr, law.rhs)


@validate_input(
    temperature_=temperature,
    entropy_=entropy,
    chemical_potential_=chemical_potential,
    particle_count_=particle_count,
)
@validate_output(enthalpy)
def calculate_enthalpy(
    temperature_: Quantity,
    entropy_: Quantity,
    chemical_potential_: Quantity,
    particle_count_: int,
) -> Quantity:
    # Note that since entropy is known up to a constant, enthalpy is also only known
    # up to a constant.

    result = law.rhs.subs({
        temperature: temperature_,
        entropy: entropy_,
        chemical_potential: chemical_potential_,
        particle_count: particle_count_,
    })
    return Quantity(result)
