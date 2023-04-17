from sympy import Expr
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Spring accumulates energy while being deformated. This law is known as Hooke's law.
## Law: E = k * x**2 / 2, where
## E is potential energy of deformated spring
## k is elastic koefficient
## x is deformation

# Conditions.
## Deformation is elactic (reversible).

spring_energy = Symbol("spring_energy", units.energy)
elastic_koefficient = Symbol("elastic_koefficient", units.force / units.length)
deformation = Symbol("deformation", units.length)

law = Eq(spring_energy, elastic_koefficient * deformation**2 / 2)

def print(expr: Expr) -> str:
    symbols = [spring_energy, elastic_koefficient, deformation]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(elastic_koefficient_=elastic_koefficient, deformation_=deformation)
@validate_output_symbol(spring_energy)
def calculate_energy(elastic_koefficient_: Quantity, deformation_: Quantity) -> Quantity:
    result_energy_expr = solve(law, spring_energy, dict=True)[0][spring_energy]        
    result_expr = result_energy_expr.subs({elastic_koefficient: elastic_koefficient_, deformation: deformation_})
    return expr_to_quantity(result_expr)
