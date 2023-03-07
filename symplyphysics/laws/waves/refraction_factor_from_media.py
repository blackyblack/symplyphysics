from symplyphysics import (
    symbols, Eq, pretty, solve, units, Quantity,
    validate_input, validate_output, expr_to_quantity, sqrt
)

# Description
## How media refracts electromagnetical waves depends on how this media transfers electrical and magnetical fields.

# Law: n = sqrt(epsilon * mu), where
## n is refraction factor of media,
## epsilon is relative dielectric permeability,
## mu is relative magnetic permeability.

refraction_factor = symbols("refraction_factor")
relative_dielectric_permeability = symbols("relative_dielectric_permeability")
relative_magnetic_permeability = symbols("relative_magnetic_permeability")

law = Eq(refraction_factor, sqrt(relative_dielectric_permeability * relative_magnetic_permeability))

def print():
    return pretty(law, use_unicode=False)

@validate_input()
@validate_output()
def calculate_refraction_factor(epsilon_: Quantity, mu_: Quantity) -> Quantity:        
    result_expr = solve(law, refraction_factor, dict=True)[0][refraction_factor]
    factor_applied = result_expr.subs({relative_dielectric_permeability: epsilon_, 
                                          relative_magnetic_permeability: mu_})
    return expr_to_quantity(factor_applied, 'refraction_factor')
