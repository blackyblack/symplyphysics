from sympy import Eq, Derivative, Point2D
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.geometry.line import two_point_function

# Description
## Gibbs-Helmholtz relations are a set of equations that relate thermodynamic potentials between each other.

# Law: U = F - T * (dF/dT)_V
## U - internal energy
## F - Helmholtz free energy
## T - absolute temperature
## V - volume

internal_energy = Symbol("internal_energy", units.energy)
free_energy = Function("free_energy", units.energy)
temperature = symbols.thermodynamics.temperature
volume = Symbol("volume", units.volume)

# Volume is held constant during the evalution of the derivative

law = Eq(
    internal_energy,
    free_energy(temperature, volume) - temperature * Derivative(free_energy(temperature, volume), temperature)
)

# TODO: Derive from definition of free energy and thermodynamical Maxwell relations


def print_law() -> str:
    return print_expression(law)


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

    result = law.rhs.subs(
        free_energy(temperature, volume), free_energy_
    ).doit().subs(
        temperature, temperature_
    )

    return Quantity(result)
