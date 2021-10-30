from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## The molar mass of a chemical compound is defined as the mass of a sample of that compound divided
## by the amount of substance in that sample, measured in moles.

## Law: M = m / mole_count
## Where:
## M (atomic or molecular weight) is the total weight of an atom or molecule
## m is the mass of the substance
## mole_count is amount of the substance, measured in moles.

atomic_weight, substance_mass, mole_count = symbols('atomic_weight substance_mass mole_count')
law = Eq(atomic_weight, substance_mass / mole_count)

def print():
    return pretty(law, use_unicode=False)

@validate_input(substance_mass_=units.mass, mole_count_=units.amount_of_substance)
@validate_output(units.mass / units.amount_of_substance)
def calculate_atomic_weight(substance_mass_: Quantity, mole_count_: Quantity) -> Quantity:
    solved = solve(law, atomic_weight, dict=True)[0][atomic_weight]
    result_expr = solved.subs({substance_mass: substance_mass_, mole_count: mole_count_})
    return expr_to_quantity(result_expr, 'atomic_weight')