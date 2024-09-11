"""
Internal energy via Helmholtz free energy
=========================================

*Gibbs—Helmholtz relations* are a set of equations that relate thermodynamic potentials between each other.

**Conditions:**

#. The number of particles in the system is held constant.
"""

from sympy import Eq, Derivative, Point2D, solve
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)
from symplyphysics.core.geometry.line import two_point_function
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import (
    helmholtz_free_energy_via_internal_energy as free_energy_def,
    entropy_is_derivative_of_free_energy as entropy_law,
)

internal_energy = Symbol("internal_energy", units.energy)
"""
Internal energy of the system.

Symbol:
    :code:`U`
"""

free_energy = Function("free_energy", units.energy)
"""
Helmholtz free energy of the system as a function of temperature and volume.

Symbol:
    :code:`F(T, V)`
"""

temperature = symbols.temperature
"""
:attr:`~symplyphysics.symbols.temperature` of the system.
"""

volume = Symbol("volume", units.volume)
"""
Volume of the system.

Symbol:
    :code:`V`
"""

# Volume is held constant during the evalution of the derivative

law = Eq(
    internal_energy,
    free_energy(temperature, volume) -
    temperature * Derivative(free_energy(temperature, volume), temperature))
r"""
:code:`U = F(T, V) - T * Derivative(F(T, V), T)`

Latex:
    .. math::
        U = F(T, V) - T \left( \frac{\partial F}{\partial T} \right)_V
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
