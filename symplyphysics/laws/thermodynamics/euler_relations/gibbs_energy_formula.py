from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import isobaric_reaction_potential as gibbs_energy_def
from symplyphysics.laws.thermodynamics.euler_relations import enthalpy_formula

# Description
## Gibbs energy differential cannot be integrated directly using the Euler's theorem on homogeneous
## functions but it can be derived via its definition and the relation for internal energy.

# Law: G = mu * N
## G - Gibbs energy
## mu - chemical potential
## N - particle count

# Note
## - This formula works for single-component system. For many-component systems replace the RHS with
##   a sum over each type of component.

gibbs_energy = Symbol("gibbs_energy", units.energy)
chemical_potential = Symbol("chemical_potential", units.energy)
particle_count = Symbol("particle_count", dimensionless)

law = Eq(gibbs_energy, chemical_potential * particle_count)

# Derive from Gibbs energy definition and enthalpy formula

_enthalpy_expr = enthalpy_formula.law.rhs.subs({
    enthalpy_formula.temperature: gibbs_energy_def.temperature,
    enthalpy_formula.entropy: gibbs_energy_def.entropy,
    enthalpy_formula.chemical_potential: chemical_potential,
    enthalpy_formula.particle_count: particle_count,
})

_gibbs_energy_expr = gibbs_energy_def.law.rhs.subs({
    gibbs_energy_def.thermal_effect: _enthalpy_expr,
})

assert expr_equals(_gibbs_energy_expr, law.rhs)


def print_law() -> str:
    return print_expression(law)


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
