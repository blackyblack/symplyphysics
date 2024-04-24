from sympy import Eq, Derivative, Point2D
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

# Description
## Gibbs-Helmholtz relations are a set of equations that relate thermodynamic potentials between each other.

# Law: H = G - T * (dG/dT)_p
## H - enthalpy
## G - Gibbs energy
## T - absolute temperature
## p - pressure
## (d/dT)_p - derivative with respect to temperature at constant pressure

enthalpy = Symbol("enthalpy", units.energy)
gibbs_energy = Function("gibbs_energy", units.energy)
temperature = symbols.thermodynamics.temperature
pressure = Symbol("pressure", units.pressure)

law = Eq(
    enthalpy,
    gibbs_energy(temperature, pressure) - temperature * Derivative(gibbs_energy(temperature, pressure), temperature)
)

# TODO: Derive from definition of Gibbs energy and thermodynamical relations


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
        gibbs_energy(temperature, pressure), gibbs_energy_,
    ).doit().subs(
        temperature, temperature_
    )

    return Quantity(result)
