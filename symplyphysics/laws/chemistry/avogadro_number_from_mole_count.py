from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units, S,
    validate_input, expr_to_quantity, convert_to
)

# Description
## The Avogadro constant is the proportionality factor that relates the number of constituent particles
## (usually molecules, atoms or ions) in a sample with the amount of substance in that sample.

## Law: Na = n / mole_count
## Where:
## n is the number of particles
## mole_count is the number of moles of the substance
## Na is Avogadro's number

particles_count, mole_count = symbols('particles_count mole_count')
law = Eq(units.avogadro, particles_count / mole_count)

def print():
    return pretty(law, use_unicode=False)

@validate_input(mole_count_=units.amount_of_substance)
def calculate_particles_count(mole_count_: Quantity) -> int:
    solved = solve(law, particles_count, dict=True)[0][particles_count]
    result_expr = solved.subs({mole_count: mole_count_})
    result_quant = expr_to_quantity(result_expr, 'particles_count')
    return int(convert_to(result_quant, S.One).n())