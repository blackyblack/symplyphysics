from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols)

# Description
## The molar mass of a chemical compound is defined as the mass of a sample of that compound divided
## by the amount of substance in that sample, measured in moles.

## Law: M = m / mole_count
## Where:
## M (atomic or molecular weight) is the total weight of an atom or molecule
## m is the mass of the substance
## mole_count is amount of the substance, measured in moles.

atomic_weight = Symbol("atomic_weight", units.mass / units.amount_of_substance)
mole_count = Symbol("mole_count", units.amount_of_substance)

law = Eq(atomic_weight, symbols.basic.mass / mole_count)


def print_law() -> str:
    return print_expression(law)


@validate_input(substance_mass_=symbols.basic.mass, mole_count_=mole_count)
@validate_output(atomic_weight)
def calculate_atomic_weight(substance_mass_: Quantity, mole_count_: Quantity) -> Quantity:
    solved = solve(law, atomic_weight, dict=True)[0][atomic_weight]
    result_expr = solved.subs({symbols.basic.mass: substance_mass_, mole_count: mole_count_})
    return Quantity(result_expr)
