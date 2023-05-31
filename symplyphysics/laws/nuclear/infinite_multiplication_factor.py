from sympy import (Eq, solve, symbols)
from symplyphysics import print_expression
from symplyphysics.core.symbols.probability import Probability

# Description
## Infinite multiplication factor: k_infinite = η * ε * p * f
## Where:
## η - neutron reproduction factor.
##   See [neutron reproduction factor](./reproduction_factor_from_macroscopic_fission_cross_section.py) implementation.
## ε - fast fission factor.
##   See [fast fission factor](./fast_fission_factor_from_resonance_escape_probability.py) implementation.
## p - resonance escape probability.
##   See [resonance escape probability](./resonance_escape_probability_from_resonance_absorption_integral.py) implementation.
## f - thermal utilisation factor.
##   See [thermal utilisation factor](./thermal_utilisation_factor_from_macroscopic_absorption_cross_sections.py) implementation.
## k_infinite (infinite multiplication factor) is the ratio of the neutrons produced by fission in one neutron
##   generation to the number of neutrons lost through absorption in the preceding neutron generation

neutron_reproduction = symbols("neutron_reproduction")
thermal_utilisation = symbols("thermal_utilisation")
resonance_escape_probability = symbols("resonance_escape_probability")
fast_fission = symbols("fast_fission")
infinite_multiplication_factor = symbols("infinite_multiplication_factor")

law = Eq(infinite_multiplication_factor,
    neutron_reproduction * fast_fission * resonance_escape_probability * thermal_utilisation)


def print() -> str:
    return print_expression(law)


def calculate_multiplication_factor(neutron_reproduction_: float, fast_fission_: float,
    resonance_escape_probability_: Probability, thermal_utilisation_: Probability) -> float:

    result_factor_expr = solve(law, infinite_multiplication_factor,
        dict=True)[0][infinite_multiplication_factor]
    result_expr = result_factor_expr.subs({
        neutron_reproduction: neutron_reproduction_,
        fast_fission: fast_fission_,
        resonance_escape_probability: resonance_escape_probability_.value,
        thermal_utilisation: thermal_utilisation_.value
    })
    return result_expr.evalf()
