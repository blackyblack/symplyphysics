from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Enthalpy of a thermodynamic system is defined as the sum of its internal energy and the product of
## its pressure and volume, which is sometimes referred to as the pressure energy.

# Law: H = U + p * V
## H - enthalpy of the system
## U - internal energy of the system
## p - pressure inside the system
## V - volume of the system

enthalpy = Symbol("enthalpy", units.energy)
internal_energy = Symbol("internal_energy", units.energy)
pressure = Symbol("pressure", units.pressure)
volume = Symbol("volume", units.volume)

law = Eq(enthalpy, internal_energy + pressure * volume)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    internal_energy_=internal_energy,
    pressure_=pressure,
    volume_=volume,
)
@validate_output(enthalpy)
def calculate_enthalpy(
    internal_energy_: Quantity,
    pressure_: Quantity,
    volume_: Quantity,
) -> Quantity:
    # Note that technically the internal energy is only known up to a constant

    result = law.rhs.subs({
        internal_energy: internal_energy_,
        pressure: pressure_,
        volume: volume_,
    })
    return Quantity(result)
