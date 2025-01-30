"""
Intensive parameters relation
=============================

The *Gibbs—Duhem relation* is a relationship among the intensive parameters of the system.
Subsequently, for a system with :math:`i` components, there are :math:`(i + 1)` independent
parameters, or degrees of freedom.

**Notation:**

#. :math:`d` denotes an exact, path-independent differential.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Gibbs%E2%80%93Duhem_equation>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import internal_energy_differential
from symplyphysics.laws.thermodynamics.euler_relations import internal_energy_formula

entropy = symbols.entropy
"""
:symbols:`entropy` of the system.
"""

temperature_change = clone_as_symbol(symbols.temperature, display_symbol="dT", display_latex="dT")
"""
Infinitesimal change in :symbols:`temperature` of the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

pressure_change = clone_as_symbol(symbols.pressure, display_symbol="dp", display_latex="dp")
"""
Infinitesimal change in :symbols:`pressure` inside the system.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

chemical_potential_change = clone_as_symbol(symbols.chemical_potential, display_symbol="d(mu)", display_latex="d\\mu")
"""
Infinitesimal change in :symbols:`chemical_potential` of the system.
"""

law = Eq(
    entropy * temperature_change - volume * pressure_change +
    particle_count * chemical_potential_change,
    0,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from internal energy representations

_internal_energy_diff_eqn = internal_energy_differential.law.subs({
    internal_energy_differential.temperature: internal_energy_formula.temperature,
    internal_energy_differential.pressure: internal_energy_formula.pressure,
    internal_energy_differential.chemical_potential: internal_energy_formula.chemical_potential,
})

_internal_energy_expr = internal_energy_formula.law.rhs.subs({
    internal_energy_formula.entropy: entropy,
    internal_energy_formula.volume: volume,
    internal_energy_formula.particle_count: particle_count,
})

_internal_energy_diff_expr = (
    _internal_energy_expr.diff(internal_energy_formula.temperature) * temperature_change +
    _internal_energy_expr.diff(entropy) * internal_energy_differential.entropy_change +
    _internal_energy_expr.diff(internal_energy_formula.pressure) * pressure_change +
    _internal_energy_expr.diff(volume) * internal_energy_differential.volume_change +
    _internal_energy_expr.diff(internal_energy_formula.chemical_potential) *
    chemical_potential_change +
    _internal_energy_expr.diff(particle_count) * internal_energy_differential.particle_count_change)

_chemical_potential_change_derived = solve(
    (_internal_energy_diff_eqn,
    Eq(internal_energy_differential.internal_energy_change, _internal_energy_diff_expr)),
    (internal_energy_differential.internal_energy_change, chemical_potential_change),
    dict=True,
)[0][chemical_potential_change]

_chemical_potential_change_from_law = solve(law, chemical_potential_change)[0]

assert expr_equals(_chemical_potential_change_derived, _chemical_potential_change_from_law)


@validate_input(
    entropy_=entropy,
    temperature_change_=temperature_change,
    volume_=volume,
    pressure_change_=pressure_change,
    particle_count_=particle_count,
)
@validate_output(chemical_potential_change)
def calculate_chemical_potential_change(
    entropy_: Quantity,
    temperature_change_: Quantity,
    volume_: Quantity,
    pressure_change_: Quantity,
    particle_count_: int,
) -> Quantity:
    result = _chemical_potential_change_from_law.subs({
        entropy: entropy_,
        temperature_change: temperature_change_,
        volume: volume_,
        pressure_change: pressure_change_,
        particle_count: particle_count_,
    })
    return Quantity(result)
