from symplyphysics import (
    symbols, Eq, pretty, solve, Probability
)

# Description
## Infinite multiplication factor: k_infinite = η * ε * p * f
## Where:
## η - neutron reproduction factor.
##   See [neutron reproduction factor](./reproduction_factor_from_macroscopic_fission_cross_section.py) implementation.
## ε (fast fission factor) is the ratio of the fast neutrons produced by fissions at all energies
##   to the number of fast neutrons produced in thermal fission
## p (resonance escape probability) is the ratio of the number of neutrons that reach thermal energies
##   to the number of fast neutrons that start to slow down
## f - thermal utilisation factor.
##   See [thermal utilisation factor](./thermal_utilisation_factor_from_macroscopic_absorption_cross_sections.py) implementation.
## k_infinite (infinite multiplication factor) is the ratio of the neutrons produced by fission in one neutron
##   generation to the number of neutrons lost through absorption in the preceding neutron generation

neutron_reproduction, thermal_utilisation = symbols('neutron_reproduction thermal_utilisation')
resonance_escape_probability, fast_fission = symbols('resonance_escape_probability fast_fission')
infinite_multiplication_factor = symbols('infinite_multiplication_factor')

law = Eq(infinite_multiplication_factor,
    neutron_reproduction * fast_fission * resonance_escape_probability * thermal_utilisation)

def print():
    return pretty(law, use_unicode=False)

def calculate_multiplication_factor(neutron_reproduction_: float, fast_fission_: float,
    resonance_escape_probability_: Probability, thermal_utilisation_: Probability) -> float:

    result_factor_expr = solve(law, infinite_multiplication_factor, dict=True)[0][infinite_multiplication_factor]
    result_expr = result_factor_expr.subs({
        neutron_reproduction: neutron_reproduction_,
        fast_fission: fast_fission_,
        resonance_escape_probability: resonance_escape_probability_.value,
        thermal_utilisation: thermal_utilisation_.value})
    return result_expr.evalf()
