from sympy import Expr
from symplyphysics import (
    assert_equivalent_dimension, Eq, pretty, solve, units, S,
    expr_to_quantity, convert_to
)
from symplyphysics.core.quantity_decorator import validate_input_symbols
from symplyphysics.core.symbols.quantities import Dimensionless, Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

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

def print(expr: Expr) -> str:
    symbols = [particles_count, mole_count]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(mole_count_=mole_count)
def calculate_particles_count(mole_count_: Quantity) -> int:
    solved = solve(law, particles_count, dict=True)[0][particles_count]
    result_expr = solved.subs(mole_count, mole_count_)
    result = expr_to_quantity(result_expr)
    assert_equivalent_dimension(result, "validate_output", "return", "calculate_particles_count", particles_count.dimension)
    return int(convert_to(result, S.One).evalf())
