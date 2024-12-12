"""
Internal energy via Helmholtz free energy
=========================================

*Gibbs—Helmholtz relations* are a set of equations that relate thermodynamic potentials between each other.

**Conditions:**

#. The number of particles in the system is held constant.

**Links:**

#. `Wikipedia, see equivalent form of this law in table <https://en.wikipedia.org/wiki/Gibbs%E2%80%93Helmholtz_equation>`__.
"""

from sympy import Eq, Derivative, Point2D, solve
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    clone_as_function,
)
from symplyphysics.core.geometry.line import two_point_function
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import (
    helmholtz_free_energy_via_internal_energy as free_energy_def,
    entropy_is_derivative_of_free_energy as entropy_law,
)

internal_energy = symbols.internal_energy
"""
:symbols:`internal_energy` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

free_energy = clone_as_function(symbols.helmholtz_free_energy, [temperature, volume])
"""
:symbols:`helmholtz_free_energy` of the system as a function of temperature and volume.
"""

# Volume is held constant during the evalution of the derivative

law = Eq(
    internal_energy,
    free_energy(temperature, volume) -
    temperature * Derivative(free_energy(temperature, volume), temperature))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from definition of free energy and thermodynamical Maxwell relations

_entropy_expr = entropy_law.law.rhs.subs({
    entropy_law.temperature: temperature,
    entropy_law.volume: volume,
}).subs(
    entropy_law.free_energy(temperature, volume, entropy_law.particle_count),
    free_energy(temperature, volume),
)

_internal_energy_expr = solve(free_energy_def.law, free_energy_def.internal_energy)[0].subs({
    free_energy_def.helmholtz_free_energy: free_energy(temperature, volume),
    free_energy_def.temperature: temperature,
    free_energy_def.entropy: _entropy_expr,
})

assert expr_equals(_internal_energy_expr, law.rhs)


@validate_input(
    free_energy_before_=free_energy,
    free_energy_after_=free_energy,
    temperature_before_=temperature,
    temperature_after_=temperature,
    temperature_=temperature,
)
@validate_output(internal_energy)
def calculate_internal_energy(
    free_energy_before_: Quantity,
    free_energy_after_: Quantity,
    temperature_before_: Quantity,
    temperature_after_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    free_energy_ = two_point_function(
        Point2D(temperature_before_, free_energy_before_),
        Point2D(temperature_after_, free_energy_after_),
        temperature,
    )

    result = law.rhs.subs(free_energy(temperature, volume),
        free_energy_).doit().subs(temperature, temperature_)

    return Quantity(result)
