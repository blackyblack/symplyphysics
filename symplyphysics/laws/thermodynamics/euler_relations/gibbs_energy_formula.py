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
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import gibbs_energy_via_enthalpy as gibbs_energy_def
from symplyphysics.laws.thermodynamics.euler_relations import enthalpy_formula

gibbs_energy = symbols.gibbs_energy
"""
:symbols:`gibbs_energy` of the system.
"""

chemical_potential = symbols.chemical_potential
"""
:symbols:`chemical_potential` of the system.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

law = Eq(gibbs_energy, chemical_potential * particle_count)
"""
:laws:symbol::

:laws:latex::
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
