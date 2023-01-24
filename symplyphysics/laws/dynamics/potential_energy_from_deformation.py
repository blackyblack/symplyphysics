from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity, SI
)

# Description
## Spring accumulates energy while being deformated. This law is known as Hooke's law.
## Law: E = k * x**2 / 2, where
## E is potential energy of deformated spring
## k is elastic koefficient
## x is deformation

# Conditions.
## Deformation is elactic (reversible).

spring_energy, elastic_koefficient, deformation = symbols('spring_energy elastic_koefficient deformation')
law = Eq(spring_energy, elastic_koefficient * deformation**2 / 2)

def print():
    return pretty(law, use_unicode=False)

@validate_input(elastic_koefficient_ = units.force / units.length, deformation_ = units.length)
@validate_output(units.energy)
def calculate_energy(elastic_koefficient_: Quantity, deformation_: Quantity) -> Quantity:
    result_energy_expr = solve(law, spring_energy, dict=True)[0][spring_energy]        
    result_expr = result_energy_expr.subs({elastic_koefficient: elastic_koefficient_, deformation: deformation_})
    return expr_to_quantity(result_expr, 'spring_energy')
