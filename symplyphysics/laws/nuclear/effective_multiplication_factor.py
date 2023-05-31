from sympy import (Eq, solve, symbols)
from symplyphysics import print_expression
from symplyphysics.core.symbols.probability import Probability

# Description
## Effective multiplication factor: k_effective = k_infinite * Pf * Pt
## Where:
## Pf - fast non-leakage factor.
##   See [fast non_leakage probability](./fast_non_leakage_probability_from_fermi_age.py)
## Pt - thermal non-leakage factor.
##   See [thermal non_leakage probability](./thermal_non_leakage_probability_from_diffusion_length.py)
## k_infinite (infinite multiplication factor) is the ratio of the neutrons produced by fission in one neutron
##   generation to the number of neutrons lost through absorption in the preceding neutron generation.
##   See [infinite multiplication factor](./infinite_multiplication_factor.py) implementation.
## k_effective (effective multiplication factor) is ratio of the neutrons produced by fission in one neutron
##   generation to the number of neutrons lost through absorption and leakage in the preceding neutron generation.

fast_non_leakage_probability = symbols("fast_non_leakage_probability")
thermal_non_leakage_probability = symbols("thermal_non_leakage_probability")
infinite_multiplication_factor = symbols("infinite_multiplication_factor")
effective_multiplication_factor = symbols("effective_multiplication_factor")

law = Eq(
    effective_multiplication_factor,
    infinite_multiplication_factor * fast_non_leakage_probability * thermal_non_leakage_probability)


def print() -> str:
    return print_expression(law)


def calculate_multiplication_factor(infinite_multiplication_factor_: float,
    fast_non_leakage_probability_: Probability,
    thermal_non_leakage_probability_: Probability) -> float:

    result_factor_expr = solve(law, effective_multiplication_factor,
        dict=True)[0][effective_multiplication_factor]
    result_expr = result_factor_expr.subs({
        infinite_multiplication_factor: infinite_multiplication_factor_,
        fast_non_leakage_probability: fast_non_leakage_probability_.value,
        thermal_non_leakage_probability: thermal_non_leakage_probability_.value
    })
    return result_expr.evalf()
