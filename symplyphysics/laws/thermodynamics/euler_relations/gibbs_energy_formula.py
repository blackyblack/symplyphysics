"""
Gibbs energy formula
====================

Gibbs energy differential cannot be integrated directly using the Euler's theorem on homogeneous
functions but it can be derived via its definition and the relation for internal energy.

**Notes:**

#. This formula words for a single-component system. For multi-component system replace the
   right-hand side with a sum over each type of components.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import gibbs_energy_via_enthalpy as gibbs_energy_def
from symplyphysics.laws.thermodynamics.euler_relations import enthalpy_formula

gibbs_energy = Symbol("gibbs_energy", units.energy)
"""
Gibbs energy of the system.

Symbol:
    :code:`G`
"""

chemical_potential = Symbol("chemical_potential", units.energy)
r"""
Chemical potential of the system.

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

law = Eq(gibbs_energy, chemical_potential * particle_count)
r"""
:code:`G = mu * N`

Latex:
    .. math::
        G = \mu N
"""

# Derive from Gibbs energy definition and enthalpy formula

_enthalpy_expr = enthalpy_formula.law.rhs.subs({
    enthalpy_formula.temperature: gibbs_energy_def.temperature,
    enthalpy_formula.entropy: gibbs_energy_def.entropy,
    enthalpy_formula.chemical_potential: chemical_potential,
    enthalpy_formula.particle_count: particle_count,
})

_gibbs_energy_expr = gibbs_energy_def.law.rhs.subs({
    gibbs_energy_def.enthalpy: _enthalpy_expr,
})

assert expr_equals(_gibbs_energy_expr, law.rhs)


@validate_input(
    chemical_potential_=chemical_potential,
    particle_count_=particle_count,
)
@validate_output(gibbs_energy)
def calculate_gibbs_energy(
    chemical_potential_: Quantity,
    particle_count_: int,
) -> Quantity:
    result = law.rhs.subs({
        chemical_potential: chemical_potential_,
        particle_count: particle_count_,
    })
    return Quantity(result)
