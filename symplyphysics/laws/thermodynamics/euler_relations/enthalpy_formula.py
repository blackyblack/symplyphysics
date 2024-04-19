from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import enthalpy_is_internal_energy_plus_pressure_energy as enthalpy_def
from symplyphysics.laws.thermodynamics.euler_relations import internal_energy_formula

# Description
## The enthalpy differential can be integrated using the Euler's theorem on homogeneous functions for
## internal energy.

# Law: H = T * S + mu * N
## H - enthalpy
## T - temperature
## S - entropy
## mu - chemical potential
## N - particle count

# Note
## - This formula works for a single-component system. For multi-component systems replace the RHS with
##   a sum over each type of component.

enthalpy = Symbol("enthalpy", units.energy)
temperature = symbols.thermodynamics.temperature
entropy = Symbol("entropy", units.energy / units.temperature)
chemical_potential = Symbol("chemical_potential", units.energy)
particle_count = Symbol("particle_count", dimensionless)

law = Eq(enthalpy, temperature * entropy + chemical_potential * particle_count)

# Derive from enthalpy definition and internal energy Euler relation

_internal_energy_expr = internal_energy_formula.law.rhs.subs({
    internal_energy_formula.temperature: temperature,
    internal_energy_formula.entropy: entropy,
    internal_energy_formula.pressure: enthalpy_def.pressure,
    internal_energy_formula.volume: enthalpy_def.volume,
    internal_energy_formula.chemical_potential: chemical_potential,
    internal_energy_formula.particle_count: particle_count,
})

_enthalpy_expr = enthalpy_def.law.rhs.subs({
    enthalpy_def.internal_energy: _internal_energy_expr
})

assert expr_equals(_enthalpy_expr, law.rhs)


def print_law() -> str:
    return print_expression(law)


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
