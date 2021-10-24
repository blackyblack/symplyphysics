from symplyphysics import (
    symbols, Eq, pretty, solve, Probability
)

# Description
## Effective multiplication factor: k_effective = k_infinite * Pf * Pt
## Where:
## Pf (fast non-leakage factor) is the ratio of the number of fast neutrons that do not leak from the reactor
##   core during the slowing down process to the number of fast neutrons produced by fissions at all energies
## Pt (thermal non-leakage factor) is the ratio of the number of thermal neutrons that do not leak from the
##   reactor core during the neutron diffusion process to the number of neutrons that reach thermal energies
## k_infinite (infinite multiplication factor) is the ratio of the neutrons produced by fission in one neutron
##   generation to the number of neutrons lost through absorption in the preceding neutron generation
## k_effective (effective multiplication factor) is ratio of the neutrons produced by fission in one neutron
##   generation to the number of neutrons lost through absorption and leakage in the preceding neutron generation

fast_non_leakage_probability = symbols('fast_non_leakage_probability')
thermal_non_leakage_probability = symbols('thermal_non_leakage_probability')
infinite_multiplication_factor = symbols('infinite_multiplication_factor')
effective_multiplication_factor = symbols('effective_multiplication_factor')

law = Eq(effective_multiplication_factor,
    infinite_multiplication_factor * fast_non_leakage_probability * thermal_non_leakage_probability)

def print():
    return pretty(law, use_unicode=False)

def calculate_multiplication_factor(infinite_multiplication_factor_: float,
    fast_non_leakage_probability_: Probability, thermal_non_leakage_probability_: Probability) -> float:

    result_factor_expr = solve(law, effective_multiplication_factor, dict=True)[0][effective_multiplication_factor]
    result_expr = result_factor_expr.subs({
        infinite_multiplication_factor: infinite_multiplication_factor_,
        fast_non_leakage_probability: fast_non_leakage_probability_.value,
        thermal_non_leakage_probability: thermal_non_leakage_probability_.value})
    return result_expr.evalf()
