"""
Enthalpy formula
================

The enthalpy differential can be integrated using the Euler's theorem on homogeneous functions for
internal energy.

**Notes:**

#. This formula words for a single-component system. For multi-component system replace the
   product of chemical potential and particle count with a sum over each type of components.
"""

from sympy import Eq
from symplyphysics import (units, dimensionless, Quantity, Symbol, validate_input, validate_output,
    symbols)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import enthalpy_is_internal_energy_plus_pressure_energy as enthalpy_def
from symplyphysics.laws.thermodynamics.euler_relations import internal_energy_formula

enthalpy = Symbol("enthalpy", units.energy)
"""
:doc:`Enthalpy <laws.thermodynamics.enthalpy_is_internal_energy_plus_pressure_energy>` of the system.

Symbol:
    :code:`H`
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

entropy = Symbol("entropy", units.energy / units.temperature)
"""
:doc:`Entropy <laws.thermodynamics.entropy_change_in_reversible_process>` of the system.

Symbol:
    :code:`S`
"""

chemical_potential = Symbol("chemical_potential", units.energy)
r"""
:doc:`Chemical potential <laws.thermodynamics.chemical_potential_is_particle_count_derivative_of_internal_energy>`
of the system.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

particle_count = Symbol("particle_count", dimensionless)
"""
Number of particles in the system.

Symbol:
    :code:`N`
"""

law = Eq(enthalpy, temperature * entropy + chemical_potential * particle_count)
r"""
:code:`H = T * S + mu * N`

Latex:
    .. math::
        H = T S + \mu N
"""

# Derive from enthalpy definition and internal energy Euler relation

_internal_energy_expr = internal_energy_formula.law.rhs.subs({
    internal_energy_formula.temperature: temperature,
    internal_energy_formula.entropy: entropy,
    internal_energy_formula.pressure: enthalpy_def.pressure,
    internal_energy_formula.volume: enthalpy_def.volume,
    internal_energy_formula.chemical_potential: chemical_potential,
    internal_energy_formula.particle_count: particle_count,
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
