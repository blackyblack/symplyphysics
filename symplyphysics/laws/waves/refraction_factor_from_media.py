from symplyphysics import (
    symbols, Eq, pretty, solve, sqrt
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

def calculate_refraction_factor(epsilon_: float, mu_: float) -> float:        
    result_expr = solve(law, refraction_factor, dict=True)[0][refraction_factor]    
    factor_applied = result_expr.subs({relative_dielectric_permeability: epsilon_, 
                                       relative_magnetic_permeability: mu_})    
    return factor_applied
