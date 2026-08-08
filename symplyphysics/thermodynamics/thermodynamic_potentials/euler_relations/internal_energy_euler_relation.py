"""
Internal energy formula
=======================

The formula of the internal energy differential can be integrated using the Euler's theorem on
homogeneous functions to get the following expression.

**Notes:**

#. This formula works for a single-component system. For a multi-component system replace the
   product of chemical potential and particle count with a sum over each type of components.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Thermodynamic_potential#Euler_relations>`__.
"""

from sympy import Eq, Derivative, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.thermodynamics.thermodynamic_potentials.differentials import internal_energy_differential
from symplyphysics.thermodynamics.thermodynamic_potentials.properties import internal_energy_is_first_order_homogeneous_function as homogeneity_property

internal_energy = symbols.internal_energy
"""
:symbols:`internal_energy` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

entropy = symbols.entropy
"""
:symbols:`entropy` of the system.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` in the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

chemical_potential = symbols.chemical_potential
"""
:symbols:`chemical_potential` of the system.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

law = Eq(internal_energy,
    temperature * entropy - pressure * volume + chemical_potential * particle_count)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from the homogeneity property of internal energy and the internal energy differential
# using Euler's homogeneous function theorem.

## Internal energy is an extensive quantity: scaling the system, i.e. multiplying its extensive
## natural variables (entropy, volume, particle count) by a factor `k`, multiplies the internal
## energy by the same factor. This is expressed by the homogeneity condition
## `U(k*S, k*V, k*N) = k*U(S, V, N)`, which is exactly what allows the internal energy
## differential to be integrated into a closed-form expression for `U`.

_factor = homogeneity_property.factor

## Differentiate both sides of the homogeneity condition with respect to the scaling factor `k`.
## The condition holds for any value of `k`, so its derivative holds for any `k` as well, and we
## are free to evaluate it at any point; `k = 1` is chosen because there the arguments of `U`
## reduce to the unscaled variables `S`, `V`, `N`. This yields Euler's theorem for first-order
## homogeneous functions:
## `S * dU/dS + V * dU/dV + N * dU/dN = U(S, V, N)`.
_euler_theorem = Eq(
    homogeneity_property.homogeneity_condition.lhs.diff(_factor),
    homogeneity_property.homogeneity_condition.rhs.diff(_factor),
).subs(_factor, 1).doit()

## The partial derivatives of internal energy with respect to its natural variables are read
## off the internal energy differential by setting the other differentials to zero.
_internal_energy_diff = internal_energy_differential.law.rhs

_energy_function = homogeneity_property.internal_energy(
    homogeneity_property.entropy,
    homogeneity_property.volume,
    homogeneity_property.particle_count,
)

_euler_theorem = _euler_theorem.subs({
    Derivative(_energy_function, homogeneity_property.entropy):
        _internal_energy_diff.subs({
        internal_energy_differential.entropy_change: 1,
        internal_energy_differential.volume_change: 0,
        internal_energy_differential.particle_count_change: 0,
        }),
    Derivative(_energy_function, homogeneity_property.volume):
        _internal_energy_diff.subs({
        internal_energy_differential.entropy_change: 0,
        internal_energy_differential.volume_change: 1,
        internal_energy_differential.particle_count_change: 0,
        }),
    Derivative(_energy_function, homogeneity_property.particle_count):
        _internal_energy_diff.subs({
        internal_energy_differential.entropy_change: 0,
        internal_energy_differential.volume_change: 0,
        internal_energy_differential.particle_count_change: 1,
        }),
})

## The right-hand side of the Euler theorem equation is `U(S, V, N)`, the internal energy of the
## unscaled system, so the equation is solved for it to obtain the closed-form expression.
_internal_energy_derived = solve(
    _euler_theorem.subs(_energy_function, internal_energy),
    internal_energy,
)[0].subs({
    internal_energy_differential.temperature: temperature,
    internal_energy_differential.pressure: pressure,
    internal_energy_differential.chemical_potential: chemical_potential,
    homogeneity_property.entropy: entropy,
    homogeneity_property.volume: volume,
    homogeneity_property.particle_count: particle_count,
})

assert expr_equals(_internal_energy_derived, law.rhs)


@validate_input(
    temperature_=temperature,
    entropy_=entropy,
    pressure_=pressure,
    volume_=volume,
    chemical_potential_=chemical_potential,
    particle_count_=particle_count,
)
@validate_output(internal_energy)
def calculate_internal_energy(
    temperature_: Quantity,
    entropy_: Quantity,
    pressure_: Quantity,
    volume_: Quantity,
    chemical_potential_: Quantity,
    particle_count_: int,
) -> Quantity:
    result = law.rhs.subs({
        temperature: temperature_,
        entropy: entropy_,
        pressure: pressure_,
        volume: volume_,
        chemical_potential: chemical_potential_,
        particle_count: particle_count_,
    })
    return Quantity(result)
