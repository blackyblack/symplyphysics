from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Symbol,
    Function,
    print_expression,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import internal_energy_differential as internal_energy_law

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

homogeneity_condition = Eq(
    internal_energy(factor * entropy, factor * volume, factor * particle_count),
    factor * internal_energy(entropy, volume, particle_count),
)

# Derive from the formula of internal energy differential

_internal_energy_expr = internal_energy_law.law.rhs

_lhs = homogeneity_condition.lhs.subs(
    internal_energy(factor * entropy, factor * volume, factor * particle_count),
    _internal_energy_expr.subs({
        internal_energy_law.entropy_change: factor * internal_energy_law.entropy_change,
        internal_energy_law.volume_change: factor * internal_energy_law.volume_change,
        internal_energy_law.particle_count_change: factor * internal_energy_law.particle_count_change,
    }),
)

_rhs = homogeneity_condition.rhs.subs(
    internal_energy(entropy, volume, particle_count),
    _internal_energy_expr,
)

assert expr_equals(_lhs, _rhs)


def print_law() -> str:
    return print_expression(homogeneity_condition)
