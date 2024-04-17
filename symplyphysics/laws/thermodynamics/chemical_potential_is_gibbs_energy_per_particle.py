from sympy import Eq, Function as SymFunction, Symbol as SymSymbol, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import (
    chemical_potential_is_particle_count_derivative_of_gibbs_energy as chemical_potential_law,
)

# Description
## Chemical potential of a thermodynamic system can be viewed as the Gibbs energy of the system
## per particle. Therefore, chemical potential is an intensive physical quantity, whereas Gibbs energy
## and particle count are extensive.

# Law: mu = G / N
## mu - chemical potential
## G - [Gibbs energy](./isobaric_reaction_potential.py)
## N - particle count

chemical_potential = Symbol("chemical_potential", units.energy)
gibbs_energy = Symbol("gibbs_energy", units.energy)
particle_count = Symbol("particle_count")

law = Eq(chemical_potential, gibbs_energy / particle_count)

# Derive law from the extensive property of chemical potential

_gibbs_energy = SymFunction("gibbs_energy")
_pressure = SymSymbol("pressure")
_temperature = SymSymbol("temperature")
_factor = SymSymbol("factor")

_gibbs_energy_eqn = Eq(gibbs_energy, _gibbs_energy(_temperature, _pressure, particle_count))

# Particle count and Gibbs energy are extensive quantities, whereas temperature and pressure
# are intensive, which means that the latter do not scale with particle count, unlike Gibbs energy.
_gibbs_energy_factored_eqn = _gibbs_energy_eqn.subs({
    particle_count: _factor * particle_count,
    gibbs_energy: _factor * gibbs_energy,
})

# `_factor` is arbitrary and can be chosen to be equal to `1 / particle_count`
_gibbs_energy_subs_eqn = _gibbs_energy_factored_eqn.subs(_factor, 1 / particle_count)

_gibbs_energy_expr = solve(_gibbs_energy_subs_eqn, gibbs_energy)[0]

_chemical_potential_eqn = chemical_potential_law.law.subs({
    chemical_potential_law.temperature: _temperature,
    chemical_potential_law.pressure: _pressure,
    chemical_potential_law.particle_count: particle_count,
    chemical_potential_law.chemical_potential: chemical_potential,
}).subs(
    chemical_potential_law.gibbs_energy(_temperature, _pressure, particle_count),
    _gibbs_energy_expr,
).doit()

_chemical_potential_derived = solve(
    (_gibbs_energy_subs_eqn, _chemical_potential_eqn),
    (_gibbs_energy(_temperature, _pressure, 1), chemical_potential),
    dict=True,
)[0][chemical_potential]

assert expr_equals(_chemical_potential_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    gibbs_energy_=gibbs_energy,
    particle_count_=particle_count,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    gibbs_energy_: Quantity,
    particle_count_: int,
) -> Quantity:
    result = law.rhs.subs({
        gibbs_energy: gibbs_energy_,
        particle_count: particle_count_,
    })
    return Quantity(result)
