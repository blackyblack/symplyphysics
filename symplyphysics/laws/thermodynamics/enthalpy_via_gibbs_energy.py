r"""
Enthalpy via Gibbs energy
=========================

Gibbs-Helmholtz relations are a set of equations that relate thermodynamic potentials between each other.
For example, enthalpy :math:`H` can be found using the Gibbs energy :math:`G` under isobaric conditions.

**Conditions:**

#. Particle count must be constant.
#. Pressure in the system must be constant.
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
    gibbs_energy_via_enthalpy as gibbs_energy_def,
    entropy_is_derivative_of_gibbs_energy as entropy_law,
)

enthalpy = symbols.enthalpy
"""
:symbols:`enthalpy` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` inside the system.
"""

gibbs_energy = clone_as_function(symbols.gibbs_energy, [temperature, pressure])
"""
Gibbs energy of the system.

Symbol:
    :code:`G(T, p)`
"""

law = Eq(
    enthalpy,
    gibbs_energy(temperature, pressure) -
    temperature * Derivative(gibbs_energy(temperature, pressure), temperature))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from definition of Gibbs energy and thermodynamical relations

_entropy_expr = entropy_law.law.rhs.subs({
    entropy_law.temperature: temperature,
    entropy_law.pressure: pressure,
}).subs(
    entropy_law.gibbs_energy(temperature, pressure, entropy_law.particle_count),
    gibbs_energy(temperature, pressure),
)

_enthalpy_expr = solve(gibbs_energy_def.law, gibbs_energy_def.enthalpy)[0].subs({
    gibbs_energy_def.gibbs_energy: gibbs_energy(temperature, pressure),
    gibbs_energy_def.temperature: temperature,
    gibbs_energy_def.entropy: _entropy_expr,
})

assert expr_equals(_enthalpy_expr, law.rhs)


@validate_input(
    gibbs_energy_before_=gibbs_energy,
    gibbs_energy_after_=gibbs_energy,
    temperature_before_=temperature,
    temperature_after_=temperature,
    temperature_=temperature,
)
@validate_output(enthalpy)
def calculate_enthalpy(
    gibbs_energy_before_: Quantity,
    gibbs_energy_after_: Quantity,
    temperature_before_: Quantity,
    temperature_after_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    gibbs_energy_ = two_point_function(
        Point2D(temperature_before_, gibbs_energy_before_),
        Point2D(temperature_after_, gibbs_energy_after_),
        temperature,
    )

    result = law.rhs.subs(
        gibbs_energy(temperature, pressure),
        gibbs_energy_,
    ).doit().subs(temperature, temperature_)

    return Quantity(result)
