from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol,
                           print_expression, validate_input_symbols,
                           validate_output_symbol)

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


def print() -> str:
    return print_expression(law)


@validate_input_symbols(elastic_koefficient_=elastic_koefficient,
                        deformation_=deformation)
@validate_output_symbol(spring_energy)
def calculate_energy(elastic_koefficient_: Quantity,
                     deformation_: Quantity) -> Quantity:
    result_energy_expr = solve(law, spring_energy, dict=True)[0][spring_energy]
    result_expr = result_energy_expr.subs({
        elastic_koefficient: elastic_koefficient_,
        deformation: deformation_
    })
    return expr_to_quantity(result_expr)
