from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Symbol,
    Function,
    print_expression,
)

# Description
## Internal energy is a first-order homogeneous function of its internal variables (entropy,
## volume, and particle count).

# Law: U(k * S, k * V, k * N) = k * U(S, V, N) forall k
## U = U(S, V, N) - internal energy as a function of its natural variables
## S - entropy
## V - volume
## N - particle count
## k - any (real) number

internal_energy = Function("internal_energy", units.energy)
entropy = Symbol("entropy", units.energy / units.temperature)
volume = Symbol("volume", units.volume)
particle_count = Symbol("particle_count", dimensionless)
factor = Symbol("factor", dimensionless)

law = Eq(
    internal_energy(factor * entropy, factor * volume, factor * particle_count),
    factor * internal_energy(entropy, volume, particle_count),
)


def print_law() -> str:
    return print_expression(law)
