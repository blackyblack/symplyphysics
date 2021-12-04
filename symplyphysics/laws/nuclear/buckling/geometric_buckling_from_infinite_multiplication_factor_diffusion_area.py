from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Geometric buckling can be derived from the diffusion equation: see [diffusion equation](../diffusion_equation_from_neutron_flux.py).
## Another geometric buckling definition depends on neutron flux: see [geometric buckling](./geometric_buckling_from_neutron_flux.py) implementation.

## Law: Bg^2 = (k_infinite / k_effective - 1) / L^2
## Where:
## k_infinite - infinite multiplication factor.
##   See [infinite multiplication factor](./infinite_multiplication_factor.py)
## k_effective - effective multiplication factor.
##   See [effective multiplication factor](./effective_multiplication_factor.py)
## L^2 - diffusion area.
##   See [diffusion area](./diffusion_area_from_diffusion_coefficient.py) implementation.
## Bg^2 - geometric buckling.
##   See [geometric buckling](./geometric_buckling_from_neutron_flux.py) implementation.

geometric_buckling_squared = symbols('geometric_buckling_squared')
infinite_multiplication_factor = symbols('infinite_multiplication_factor')
effective_multiplication_factor = symbols('effective_multiplication_factor')
diffusion_area = symbols('diffusion_area')

law = Eq(geometric_buckling_squared,
    (infinite_multiplication_factor / effective_multiplication_factor - 1) / diffusion_area)

def print():
    return pretty(law, use_unicode=False)

@validate_input(diffusion_area_=units.length**2)
@validate_output(1 / units.length**2)
def calculate_geometric_buckling_squared(
    infinite_multiplication_factor_: float,
    effective_multiplication_factor_: float,
    diffusion_area_: Quantity) -> Quantity:
    result_buckling_expr = solve(law, geometric_buckling_squared, dict=True)[0][geometric_buckling_squared]
    result_expr = result_buckling_expr.subs({
        infinite_multiplication_factor: infinite_multiplication_factor_,
        effective_multiplication_factor: effective_multiplication_factor_,
        diffusion_area: diffusion_area_})
    return expr_to_quantity(result_expr, 'geometric_buckling_squared')
