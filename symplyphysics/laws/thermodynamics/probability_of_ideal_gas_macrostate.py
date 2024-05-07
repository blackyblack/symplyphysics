from typing import Sequence
from sympy import Eq, Idx, factorial
from symplyphysics import (
    dimensionless,
    Symbol,
    convert_to_float,
    validate_input,
    validate_output,
    ProductIndexed,
    SymbolIndexed,
    global_index,
    assert_equal,
)
from symplyphysics.core.symbols.probability import Probability

# Description
## Suppose there is a vessel of volume `V` with `N` identical ideal gas particles whose movement is described
## by classical mechanics. Let us divide the volume of the vessel into `m` small enough cells with volumes `V_1`
## upto `V_m`. Let us for a moment also attach a number to each gas particle to be able to tell them apart.

## The _microstate_ of the gas is the complete description of (a) the particle count of each cell of the vessel
## and (b) the particles of which numbers appear in each cell. A particle's moving within the same cell does not
## change the microstate of the gas, whereas a particle's moving from one cell to another constitutes the change
## of the microstate.

## The _macrostate_ of the gas only requires the description of the particle count in each cell of the vessel.
## If we assume the particles to be identical, it doesn't make a difference, from a macroscopic standpoint, particles
## of which number appear in which cells. Hence the probability of the macrostate is larger than that of the microstate
## by the factor of the number of permutations of the particles among the cells, which is exactly the statistical weight
## of the macrostate.

# Law: P = G * Product(p_i**N_i, i)
## P - probability of macrostate
## G - statistical weight of the macrostate
## p_i - probability of finding one particle in i-th cell
## N_i - number of particles in i-th cell

# Conditions
## - There are no external fields acting on the system
## - All particles are identical
## - Gas is ideal
## - Probabilities must sum up to 1

macrostate_probability = Symbol("macrostate_probability", dimensionless)
statistical_weight = Symbol("statistical_weight", dimensionless)
one_particle_in_cell_probability = SymbolIndexed("one_particle_in_cell_probability", dimensionless)
particle_count_in_cell = SymbolIndexed("particle_count_in_cell", dimensionless)

law = Eq(
    macrostate_probability,
    statistical_weight
    * ProductIndexed(
        one_particle_in_cell_probability[global_index]**particle_count_in_cell[global_index],
        global_index,
    ),
)


@validate_input(
    one_particle_in_cell_probabilities_=one_particle_in_cell_probability,
    particle_count_in_cells_=particle_count_in_cell,
)
@validate_output(macrostate_probability)
def calculate_macrostate_probability(
    one_particle_in_cell_probabilities_: Sequence[Probability],
    particle_count_in_cells_: Sequence[int],
) -> Probability:
    if len(one_particle_in_cell_probabilities_) != len(particle_count_in_cells_):
        raise ValueError("The number of particle counts must equal the number of probabilities")

    try:
        assert_equal(sum(one_particle_in_cell_probabilities_), 1)
    except AssertionError as e:
        raise ValueError("The probabilities must sum up to 1") from e

    statistical_weight_ = factorial(sum(particle_count_in_cells_))
    for particle_count_ in particle_count_in_cells_:
        statistical_weight_ = statistical_weight_ // factorial(particle_count_)

    result = law.rhs.subs(statistical_weight, statistical_weight_)

    local_index = Idx("local_index", (1, len(one_particle_in_cell_probabilities_)))
    result = result.subs(global_index, local_index).doit()

    for idx, (probability_, particle_count_) in enumerate(zip(
        one_particle_in_cell_probabilities_,
        particle_count_in_cells_,
    ), 1):
        result = result.subs({
            one_particle_in_cell_probability[idx]: probability_,
            particle_count_in_cell[idx]: particle_count_,
        })

    return Probability(convert_to_float(result))
