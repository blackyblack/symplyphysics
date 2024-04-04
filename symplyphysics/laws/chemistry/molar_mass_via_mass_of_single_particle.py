from sympy import Eq, solve, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.chemistry import (
    atomic_weight_from_mass_mole_count as molar_mass_law,
    avogadro_number_from_mole_count as avogadro_law,
)

# Description
## Molar mass of a substance is the mass of 1 mole of particles that it is comprised of.

molar_mass = Symbol("molas_mass", units.mass / units.amount_of_substance)
particle_mass = Symbol("particle_mass", units.mass)

law = Eq(molar_mass, particle_mass * units.avogadro)

# Derive from another definition of molar mass

_number_of_particles = SymSymbol("number_of_particles", integer=True)

_amount_of_substance = solve(
    avogadro_law.law, avogadro_law.mole_count
)[0].subs(
    avogadro_law.particles_count, _number_of_particles
)

_total_mass = _number_of_particles * particle_mass

_molar_mass_derived = molar_mass_law.law.rhs.subs({
    symbols.basic.mass: _total_mass,
    molar_mass_law.mole_count: _amount_of_substance,
})

assert expr_equals(_molar_mass_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(particle_mass_=particle_mass)
@validate_output(molar_mass)
def calculate_molar_mass(particle_mass_: Quantity) -> Quantity:
    result = law.rhs.subs(particle_mass, particle_mass_)
    return Quantity(result)
