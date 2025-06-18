"""
Internal energy is first order homogeneous function
===================================================

Internal energy is a first-order homogeneous function of its internal variables (entropy, volume,
and particle count).

**Links:**

#. `Wikipedia, first equation <https://en.wikipedia.org/wiki/Thermodynamic_potential#Euler_relations>`__.
"""

from sympy import Eq
from symplyphysics import dimensionless, Symbol, symbols, clone_as_function
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import internal_energy_differential as _internal_energy_law

entropy = symbols.entropy
"""
:symbols:`entropy`.
"""

volume = symbols.volume
"""
:symbols:`volume`
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count`.
"""

internal_energy = clone_as_function(symbols.internal_energy, [entropy, volume, particle_count])
"""
:symbols:`internal_energy` as a function of :attr:`~entropy`, :attr:`~volume`, and
:attr:`~particle_count`.
"""

factor = Symbol("k", dimensionless, real=True)
"""
Dimensionless real-valued factor.
"""

homogeneity_condition = Eq(
    internal_energy(factor * entropy, factor * volume, factor * particle_count),
    factor * internal_energy(entropy, volume, particle_count),
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from the formula of internal energy differential

_internal_energy_expr = _internal_energy_law.law.rhs

_lhs = homogeneity_condition.lhs.subs(
    internal_energy(factor * entropy, factor * volume, factor * particle_count),
    _internal_energy_expr.subs({
    _internal_energy_law.entropy_change:
    factor * _internal_energy_law.entropy_change,
    _internal_energy_law.volume_change:
    factor * _internal_energy_law.volume_change,
    _internal_energy_law.particle_count_change:
    factor * _internal_energy_law.particle_count_change,
    }),
)

_rhs = homogeneity_condition.rhs.subs(
    internal_energy(entropy, volume, particle_count),
    _internal_energy_expr,
)

assert expr_equals(_lhs, _rhs)
