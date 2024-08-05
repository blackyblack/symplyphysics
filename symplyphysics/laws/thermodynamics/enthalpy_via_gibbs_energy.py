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
    isobaric_reaction_potential as gibbs_energy_def,
    entropy_is_derivative_of_gibbs_energy as entropy_law,
)

enthalpy = Symbol("enthalpy", units.energy)
"""
Enthalpy of the system.

Symbol:
    :code:`T`
"""

gibbs_energy = Function("gibbs_energy", units.energy)
"""
Gibbs energy of the system.

Symbol:
    :code:`G`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.basic.temperature` of the system.

Symbol:
    :code:`T`
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure inside the system.

Symbol:
    :code:`p`
"""

law = Eq(
    enthalpy,
    gibbs_energy(temperature, pressure) -
    temperature * Derivative(gibbs_energy(temperature, pressure), temperature))
r"""
:code:`H = G(T, p) - T * Derivative(G(T, p), T)`

Latex:
    .. math::
        H = G(T, p) - T \left( \frac{\partial G}{\partial T} \right)_p
"""

# Derive from definition of Gibbs energy and thermodynamical relations

_entropy_expr = entropy_law.law.rhs.subs({
    entropy_law.temperature: temperature,
    entropy_law.pressure: pressure,
}).subs(
    entropy_law.gibbs_energy(temperature, pressure, entropy_law.particle_count),
    gibbs_energy(temperature, pressure),
)

_enthalpy_expr = solve(gibbs_energy_def.law, gibbs_energy_def.thermal_effect)[0].subs({
    gibbs_energy_def.isobaric_potential: gibbs_energy(temperature, pressure),
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
