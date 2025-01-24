r"""
Probability of ideal gas macrostate
===================================

Suppose there is a vessel of volume :math:`V` with :math:`N` identical ideal gas particles whose movement is described
by classical mechanics. Let us divide the volume of the vessel into :math:`m` small enough cells with volumes :math:`V_1`
up to :math:`V_m`. Let us for a moment also attach a number to each gas particle to be able to tell them apart.

The *microstate* of the gas is the complete description of

#. the particle count of each cell of the vessel
#. the particles of which numbers appear in each cell.

A particle moving within the same cell does not change the microstate of the gas, but a particle moving
from one cell to another does.

The *macrostate* of the gas only requires the description of the particle count in each cell of the vessel.
If we assume the particles to be identical, it doesn't make a difference, from a macroscopic standpoint, particles
of which number appear in which cells. Hence the probability of the macrostate is larger than that of the microstate
by the factor of the number of permutations of the particles among the cells, which is exactly the statistical weight
of the macrostate.

**Conditions:**

#. There are no external fields acting on the system.
#. All particles are identical.
#. The gas is ideal.
#. Probabilities must sum up to :math:`1`.

**Links:**

#. Formula 80.8 on p. 298 of "General Course of Physics" (Obschiy kurs fiziki), vol. 1 by Sivukhin D.V. (1979).

..
    TODO find English link
"""

from typing import Sequence
from sympy import Eq, Idx, factorial
from symplyphysics import (
    dimensionless,
    convert_to_float,
    validate_output,
    ProductIndexed,
    SymbolIndexedNew,
    global_index,
    assert_equal,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension
from symplyphysics.core.symbols.probability import Probability
from symplyphysics.core.symbols.symbols import clone_as_indexed

macrostate_probability = clone_as_symbol(symbols.probability, subscript="\\text{macro}")
"""
:symbols:`probability` of the macrostate.
"""

statistical_weight = symbols.statistical_weight
"""
:symbols:`statistical_weight` of the macrostate.
"""

particle_in_cell_probability = clone_as_indexed(symbols.probability)
"""
:symbols:`probability` of finding at least one particle in cell :math:`i`.
"""

particle_count_in_cell = SymbolIndexedNew("N", dimension=dimensionless)
"""
Number of particles in cell :math:`i`.
"""

# TODO: Create law for probability of microstate and move the product there

law = Eq(
    macrostate_probability,
    statistical_weight * ProductIndexed(
    particle_in_cell_probability[global_index]**particle_count_in_cell[global_index],
    global_index,
    ),
)
r"""
..
    The printers do not work yey with `ProductIndexed`.

:code:`P_macro = G * Product(p_i^N_i, i)`

Latex:
    .. math::
        P_\text{macro} = G \prod_i p_i^{N_i}
"""


@validate_output(macrostate_probability)
def calculate_macrostate_probability(
    probabilities_and_particle_counts_: Sequence[tuple[Probability, int]],) -> Probability:
    for probability_, particle_count_ in probabilities_and_particle_counts_:
        assert_equivalent_dimension(
            probability_,
            "probability",
            "calculate_macrostate_probability",
            particle_in_cell_probability.dimension,
        )
        assert_equivalent_dimension(
            particle_count_,
            "particle_count",
            "calculate_macrostate_probability",
            particle_count_in_cell.dimension,
        )

    probabilities_, particle_counts_ = zip(*probabilities_and_particle_counts_)

    try:
        assert_equal(sum(probabilities_), 1)
    except AssertionError as e:
        raise ValueError("The probabilities must sum up to 1") from e

    # See [statistical weight law](./maxwell_boltzmann_statistics/statistical_weight_of_macrostate.py) for details
    statistical_weight_ = factorial(sum(particle_counts_))
    for particle_count_ in particle_counts_:
        statistical_weight_ /= factorial(particle_count_)

    result = law.rhs.subs(statistical_weight, statistical_weight_)

    local_index = Idx("local_index", (1, len(probabilities_)))
    result = result.subs(global_index, local_index).doit()

    for idx, (probability_, particle_count_) in enumerate(probabilities_and_particle_counts_, 1):
        result = result.subs({
            particle_in_cell_probability[idx]: probability_,
            particle_count_in_cell[idx]: particle_count_,
        })

    return Probability(convert_to_float(result))
