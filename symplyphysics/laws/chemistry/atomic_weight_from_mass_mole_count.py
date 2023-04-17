from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## The molar mass of a chemical compound is defined as the mass of a sample of that compound divided
## by the amount of substance in that sample, measured in moles.

## Law: M = m / mole_count
## Where:
## M (atomic or molecular weight) is the total weight of an atom or molecule
## m is the mass of the substance
## mole_count is amount of the substance, measured in moles.

atomic_weight = Symbol("atomic_weight", units.mass / units.amount_of_substance)
substance_mass = Symbol("substance_mass", units.mass)
mole_count = Symbol("mole_count", units.amount_of_substance)

law = Eq(atomic_weight, substance_mass / mole_count)

def print(expr: Expr) -> str:
    symbols = [atomic_weight, substance_mass, mole_count]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(substance_mass_=substance_mass, mole_count_=mole_count)
@validate_output_symbol(atomic_weight)
def calculate_atomic_weight(substance_mass_: Quantity, mole_count_: Quantity) -> Quantity:
    solved = solve(law, atomic_weight, dict=True)[0][atomic_weight]
    result_expr = solved.subs({substance_mass: substance_mass_, mole_count: mole_count_})
    return expr_to_quantity(result_expr)
