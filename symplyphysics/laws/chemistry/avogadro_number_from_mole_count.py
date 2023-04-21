from sympy import (Eq, solve, S)
from symplyphysics import (
    units, expr_to_quantity, Quantity, Symbol, print_expression, Dimensionless, convert_to,
    validate_input_symbols,
)
from symplyphysics.core.quantity_decorator import assert_equivalent_dimension

# Description
## The Avogadro constant is the proportionality factor that relates the number of constituent particles
## (usually molecules, atoms or ions) in a sample with the amount of substance in that sample.

## Law: Na = n / mole_count
## Where:
## n is the number of particles
## mole_count is the number of moles of the substance
## Na is Avogadro's number

particles_count = Symbol("particles_count", Dimensionless)
mole_count = Symbol("mole_count", units.amount_of_substance)

law = Eq(units.avogadro, particles_count / mole_count)

def print() -> str:
    return print_expression(law)

@validate_input_symbols(mole_count_=mole_count)
def calculate_particles_count(mole_count_: Quantity) -> int:
    solved = solve(law, particles_count, dict=True)[0][particles_count]
    result_expr = solved.subs(mole_count, mole_count_)
    result = expr_to_quantity(result_expr)
    assert_equivalent_dimension(result, "validate_output", "return", "calculate_particles_count", particles_count.dimension)
    return int(convert_to(result, S.One).evalf())
